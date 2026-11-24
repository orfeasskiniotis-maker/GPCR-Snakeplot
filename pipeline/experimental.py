"""Average helix heights and phases over experimental, membrane-oriented structures.

For one GPCR class:
1. Take every structure GPCRdb lists for the class in the inactive state
   (any species) and keep the best-resolution one per receptor, so a receptor
   with fifty structures (beta2AR) counts as much as one with a single entry.
   Mixing states would blur TM6, which moves outwards by ~10 A on activation.
2. Download OPM's copy of each entry. OPM has rotated and shifted it so the
   membrane normal is the z axis and z = 0 is the middle of the bilayer.
3. Match the structure's residues to the receptor's sequence (residue numbers
   usually follow UniProt; otherwise the best constant offset is used), which
   gives each C-alpha its GPCRdb generic number.
4. Measure Z_h and phi_h for each helix (build_template.helix_params) and
   average: arithmetic mean for Z_h, circular mean for phi_h.

Usage:  python experimental.py A      (prints the per-structure table)
"""
import sys

import gemmi
import numpy as np

import sources
from build_template import helix_params, tm_residues

STATE = "Inactive"
MAX_TRIES_PER_RECEPTOR = 10  # if OPM lacks the best entry, try the next best


def class_structures(gpcr_class):
    """Best-resolution inactive structures first, grouped by receptor."""
    rows = [s for s in sources.gpcrdb_structures()
            if s["class"].startswith(f"Class {gpcr_class} ") and s["state"] == STATE]
    by_receptor = {}
    for s in sorted(rows, key=lambda s: s["resolution"] or 99):
        by_receptor.setdefault(s["protein"], []).append(s)
    return by_receptor


def half_thickness(path):
    for line in open(path):
        if "1/2 of bilayer thickness" in line:
            return float(line.split(":")[1])
    return None


def chain_ca(chain):
    """{author residue number: (one-letter code, C-alpha coords)} for standard residues."""
    out = {}
    for r in chain:
        if r.seqid.icode != " " or not r.find_atom("CA", "*"):
            continue
        info = gemmi.find_tabulated_residue(r.name)
        if info is None or not info.is_amino_acid():
            continue
        out[r.seqid.num] = (info.one_letter_code.upper(), np.array(r["CA"][0].pos.tolist()))
    return out


def match_to_sequence(cas, seq):
    """Best constant offset between structure numbering and the receptor
    sequence: returns (offset, fraction of residues whose amino acid matches).
    Residues outside the receptor sequence (fusion proteins such as T4
    lysozyme, numbered 1000+) do not count against the match."""
    nums = sorted(cas)
    best = (0, 0, 0.0)   # offset, hits, fraction
    for off in range(-60, 61):
        inside = [n for n in nums if n + off in seq]
        hits = sum(1 for n in inside if seq[n + off] == cas[n][0])
        if hits > best[1] or (hits == best[1] and off == 0):
            best = (off, hits, hits / max(len(inside), 1))
    return best[0], best[2]


def measure(structure):
    """Helix parameters of one structure, or (None, reason)."""
    path = sources.opm_structure(structure["pdb_code"])
    if path is None:
        return None, "not in OPM"
    residues = sources.gpcrdb_residues(structure["protein"])
    seq = {r["sequence_number"]: r["amino_acid"] for r in residues}
    st = gemmi.read_structure(str(path))
    best = None
    for chain in st[0]:
        cas = chain_ca(chain)
        if len(cas) < 150:
            continue
        off, frac = match_to_sequence(cas, seq)
        if best is None or frac * len(cas) > best[1] * len(best[2]):
            best = (off, frac, cas)
    if best is None:
        return None, "no receptor chain"
    off, frac, cas = best
    if frac < 0.85:
        return None, f"sequence match only {frac:.0%}"
    ca = {n + off: xyz for n, (aa, xyz) in cas.items()}
    tm = tm_residues(residues, ca)
    if len(tm) < 7:
        return None, f"only {len(tm)} TM helices resolved"
    # Which way is out? TM1 starts outside the cell, so its first residues lie
    # on the outer side.
    tm1 = sorted(tm[1])
    sign = 1.0 if ca[tm1[0][1]][2] > ca[tm1[-1][1]][2] else -1.0
    params = helix_params(tm, ca, np.array([0.0, 0.0, sign]), np.zeros(3))
    for h, p in params.items():
        expected = 1 if h % 2 == 0 else -1
        if p["dir"] != expected:
            return None, f"TM{h} points the wrong way"
        if abs(p["Z50_A"]) > 25:   # R3x50 can sit ~21 A down, at the inner surface
            return None, f"TM{h} x50 at {p['Z50_A']:.0f} A from the membrane centre"
    return {"helices": params, "half": half_thickness(path), "offset": off, "match": frac}, None


def circ_mean_sd(degrees):
    a = np.radians(degrees)
    s, c = np.mean(np.sin(a)), np.mean(np.cos(a))
    r = min(1.0, np.hypot(s, c))
    return float(np.degrees(np.arctan2(s, c))), float(np.degrees(np.sqrt(-2 * np.log(max(r, 1e-9)))))


def average_params(gpcr_class, verbose=True):
    used, skipped = [], []
    for receptor, entries in class_structures(gpcr_class).items():
        for s in entries[:MAX_TRIES_PER_RECEPTOR]:
            try:
                result, why = measure(s)
            except Exception as e:   # one bad file must not stop the class
                result, why = None, f"{type(e).__name__}: {e}"
            if result:
                used.append((s, result))
                if verbose:
                    print(f"  {s['pdb_code']} {receptor:14s} {s['resolution']} A  offset {result['offset']:+d}")
                break
            skipped.append((s["pdb_code"], receptor, why))
    helices = {}
    for h in range(1, 8):
        vals = [r["helices"][h] for _, r in used if h in r["helices"]]
        if not vals:
            continue
        z = np.array([v["Z50_A"] for v in vals])
        phi, phi_sd = circ_mean_sd([v["phi"] for v in vals])
        helices[h] = {"Z50_A": float(z.mean()), "Z50_sd": float(z.std()), "phi": phi, "phi_sd": phi_sd,
                      "dir": vals[0]["dir"], "n": len(vals)}
    halves = [r["half"] for _, r in used if r["half"]]
    if verbose:
        print(f"class {gpcr_class}: {len(used)} structures used, {len(skipped)} skipped")
        for code, rec, why in skipped:
            print(f"  skipped {code} {rec}: {why}")
    return {
        "n": len(used) if len(helices) == 7 else 0,
        "helices": helices,
        "half_thickness_A": float(np.mean(halves)) if halves else 15.0,
        "structures": [{"pdb": s["pdb_code"], "receptor": s["protein"], "resolution": s["resolution"]}
                       for s, _ in used],
    }


if __name__ == "__main__":
    res = average_params(sys.argv[1] if len(sys.argv) > 1 else "A")
    for h, p in res["helices"].items():
        print(f"TM{h}: Z50 {p['Z50_A']:+.1f} +- {p['Z50_sd']:.1f} A, phi {p['phi']:+.0f} +- {p['phi_sd']:.0f} deg (n {p['n']})")
    print(f"membrane half-thickness {res['half_thickness_A']:.1f} A")
