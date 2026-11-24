"""Build a class template: generic number -> fixed (x, y) slot on the snake plot.

For a residue on helix h at position p (see the roadmap):

    theta = phi_h + 100 deg * d_h * (p - 50)        (angle around the helix)
    x     = X_h + R * sin(theta)                    (column + helix wiggle)
    z     = Z_h + 1.5 A * d_h * (p - 50)            (height in the membrane)

d_h is +1 when the helix runs from the inside of the cell to the outside
(TM2, 4, 6) and -1 when it runs inwards (TM1, 3, 5, 7). Because the helix
direction flips, its residues also turn the other way around the membrane
normal, which is why d_h appears in theta too (the roadmap formula leaves
this sign out).

theta = 0 means the residue faces the centre of the helix bundle (towards
the ligand pocket); theta = 180 faces the lipids. The page draws residues
with cos(theta) > 0 as "front" and the rest slightly faded.

Where Z_h and phi_h come from
-----------------------------
Two sources:

experimental (used when a class has enough structures)
    Every inactive-state structure of the class that GPCRdb lists and OPM
    has placed in a membrane (best resolution per receptor, any species).
    OPM puts the membrane normal on the z axis with z = 0 at the bilayer
    centre, so for each structure:
      Z_h   = height of the x50 C-alpha (positive = outside the cell)
      phi_h = angle of residues x40-x60 around their helix axis, relative
              to the bundle centre (circular mean)
    The template takes the mean Z_h and the circular mean phi_h over all
    structures, and OPM's mean membrane half-thickness. See experimental.py.

alphafold (fallback when a class has too few structures)
    One AlphaFold model, with the membrane estimated from the model itself:
      membrane normal = average direction of the seven TM helix axes
      membrane centre = average position of the seven x50 C-alpha atoms

Usage:
    python build_template.py adrb2_human A                 # experimental, falls back to AlphaFold
    python build_template.py t2r46_human T2 --alphafold    # AlphaFold only
"""
import json
import sys
from pathlib import Path

import gemmi
import numpy as np

import sources
from generic_numbers import parse

OUT = Path(__file__).parent.parent / "site" / "data" / "templates"

# Drawing constants (SVG pixels). y = 0 is the membrane centre, up is outside.
PX_PER_A = 12.0        # vertical scale: 1.5 A rise per residue -> 18 px
RISE_A = 1.5           # alpha-helix rise per residue
TURN_DEG = 100.0       # alpha-helix turn per residue
HELIX_R = 20.0         # drawn helix radius
COLUMN = 82.0          # distance between helix columns
MEMBRANE_HALF_A = 15.0 # half the hydrophobic thickness when there is no OPM value (~30 A)
MIN_STRUCTURES = 4     # fewer usable experimental structures than this -> AlphaFold


def ca_coords(cif_path):
    st = gemmi.read_structure(str(cif_path))
    chain = st[0][0]
    return {r.seqid.num: np.array(r["CA"][0].pos.tolist()) for r in chain}


def helix_axis(points):
    """Best-fit line through C-alpha atoms: centre point and unit direction."""
    centre = points.mean(axis=0)
    direction = np.linalg.svd(points - centre)[2][0]
    if np.dot(points[-1] - points[0], direction) < 0:  # point along the sequence
        direction = -direction
    return centre, direction


def tm_residues(residues, ca):
    """{helix: [(position, sequence_number), ...]} for TM residues that have a C-alpha."""
    tm = {}
    for r in residues:
        g = parse(r["display_generic_number"])
        if r["protein_segment"].startswith("TM") and g and r["sequence_number"] in ca:
            tm.setdefault(g.helix, []).append((g.position, r["sequence_number"]))
    return tm


def helix_params(tm, ca, normal, centre):
    """Measure Z_h and phi_h of each helix.

    normal: unit membrane normal pointing out of the cell; centre: a point on
    the membrane mid-plane. Returns {h: {"Z50_A", "phi", "dir"}} for every
    helix with at least 8 residues and its x50 present.
    """
    axes = {h: helix_axis(np.array([ca[s] for _, s in v])) for h, v in tm.items() if len(v) >= 8}
    bundle = np.mean([axes[h][0] for h in axes], axis=0)
    out = {}
    for h, v in sorted(tm.items()):
        x50 = next((s for p, s in v if p == 50), None)
        if h not in axes or x50 is None:
            continue
        c, u = axes[h]
        d = 1 if np.dot(u, normal) > 0 else -1

        def angle(pos):
            """Angle (deg) of a C-alpha around the helix axis, 0 = towards bundle centre."""
            a = pos - c; a -= np.dot(a, u) * u
            t = bundle - c; t -= np.dot(t, u) * u
            return np.degrees(np.arctan2(np.dot(np.cross(t, a), normal), np.dot(t, a)))

        # Remove the ideal 100 deg/residue turn, then average what is left.
        offsets = [np.radians(angle(ca[s]) - TURN_DEG * d * (p - 50))
                   for p, s in v if 40 <= p <= 60 and p == int(p)]
        out[h] = {
            "Z50_A": float(np.dot(ca[x50] - centre, normal)),
            "phi": float(np.degrees(np.arctan2(np.mean(np.sin(offsets)), np.mean(np.cos(offsets))))),
            "dir": d,
        }
    return out


def alphafold_params(entry_name):
    """Helix parameters measured on one AlphaFold model (membrane estimated from the model)."""
    residues = sources.gpcrdb_residues(entry_name)
    accession = sources.gpcrdb_protein(entry_name)["accession"]
    model_path, _ = sources.alphafold_model(accession)
    ca = ca_coords(model_path)
    tm = tm_residues(residues, ca)
    axes = {h: helix_axis(np.array([ca[s] for _, s in v])) for h, v in tm.items()}
    # TM1 runs inwards, TM2 outwards, ... so flip every other axis and average.
    normal = sum((-1) ** h * axes[h][1] for h in axes)
    normal /= np.linalg.norm(normal)
    centre = np.mean([ca[next(s for p, s in v if p == 50)] for v in tm.values()], axis=0)
    return helix_params(tm, ca, normal, centre), f"AlphaFold model of {accession}"


def make_slots(helices, membrane_half_a):
    # Slots for every position x5-x95 (GPCRdb helices can run that far), plus a bulge slot after each one
    # (x55 -> bulge 551 at 55.1), so any receptor's bulges have a place.
    slots = {}
    for h, prm in helices.items():
        # (bulges only after two-digit positions: "x51" must stay position 51, not 5.1)
        for p in list(range(5, 96)) + [q + 0.1 for q in range(10, 95)]:
            theta = prm["phi"] + TURN_DEG * prm["dir"] * (p - 50)
            z = prm["Z50_A"] + RISE_A * prm["dir"] * (p - 50)
            key = f"{h}x{int(p)}" if p == int(p) else f"{h}x{int(p)}{round((p % 1) * 10)}"
            slots[key] = {
                "x": round(prm["X"] + HELIX_R * np.sin(np.radians(theta)), 1),
                "y": round(-z * PX_PER_A, 1),
                "helix": h,
                "front": bool(np.cos(np.radians(theta)) > 0),
            }

    # Helix 8 lies flat along the inner membrane surface, to the right of TM7.
    h8_y = (membrane_half_a + 4) * PX_PER_A
    for p in range(40, 86):   # Class B1 helix 8 runs to about 8x72
        theta = TURN_DEG * (p - 47)
        slots[f"8x{p}"] = {
            "x": round(helices[7]["X"] + 60 + (p - 47) * RISE_A * PX_PER_A, 1),
            "y": round(h8_y + 0.5 * HELIX_R * np.sin(np.radians(theta)), 1),
            "helix": 8,
            "front": bool(np.cos(np.radians(theta)) > 0),
        }
    return slots


def build(entry_name, gpcr_class, force_alphafold=False):
    measured, half_a, structures = None, MEMBRANE_HALF_A, []
    if not force_alphafold:
        import experimental
        result = experimental.average_params(gpcr_class)
        if result["n"] >= MIN_STRUCTURES:
            measured, half_a, structures = result["helices"], result["half_thickness_A"], result["structures"]
            method = "experimental"
            source = (f"Mean over {result['n']} inactive-state experimental structures "
                      f"oriented in the membrane by OPM (best resolution per receptor)")
        else:
            print(f"only {result['n']} usable inactive structures for class {gpcr_class}: using AlphaFold")
    if measured is None:
        measured, source = alphafold_params(entry_name)
        method = "alphafold"
        source += "; too few inactive experimental structures for an average"

    helices = {}
    for h in range(1, 8):
        m = measured[h]
        helices[h] = {"X": (h - 1) * COLUMN, "Z50_A": round(m["Z50_A"], 2), "phi": round(m["phi"], 1),
                      "dir": m["dir"], **{k: m[k] for k in ("n", "Z50_sd", "phi_sd") if k in m}}

    template = {
        "class": gpcr_class,
        "reference": entry_name,
        "method": method,
        "source": source,
        "structures": structures,
        "units": "SVG px; y = 0 is the membrane centre, negative y is outside the cell",
        "px_per_A": PX_PER_A,
        "membrane_half_A": round(half_a, 1),
        "membrane_half_px": round(half_a * PX_PER_A, 1),
        "helices": helices,
        "slots": make_slots(helices, half_a),
    }
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"class{gpcr_class}.json"
    out.write_text(json.dumps(template, indent=1))
    print(f"wrote {out} ({method}, {len(template['slots'])} slots, membrane half-thickness {half_a:.1f} A)")
    for h, prm in helices.items():
        spread = f"   n {prm['n']}, sd {prm['Z50_sd']:.1f} A / {prm['phi_sd']:.0f} deg" if "n" in prm else ""
        print(f"  TM{h}: Z50 {prm['Z50_A']:+.1f} A, phase {prm['phi']:+.0f} deg, dir {prm['dir']:+d}{spread}")


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    build(args[0] if args else "adrb2_human", args[1] if len(args) > 1 else "A",
          force_alphafold="--alphafold" in sys.argv)
