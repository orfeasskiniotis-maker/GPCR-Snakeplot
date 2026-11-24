"""Build the data the page needs for one receptor.

Writes:
  site/data/receptors/<entry_name>.json   residues, generic numbers, pLDDT, metadata
  site/data/models/<AF id>.cif            the AlphaFold model, served next to the page

Usage:  python build_receptor.py adrb2_human
"""
import html
import json
import re
import sys
from pathlib import Path

import gemmi
import numpy as np

import sources
from generic_numbers import parse
from build_template import ca_coords, membrane_frame
from superpose import superpose

# The receptor each class template was built from (and models are superposed on).
REFERENCE = {
    "A": "adrb2_human",    # β2-adrenoceptor
    "B1": "glp1r_human",   # GLP-1 receptor
    "B2": "agrl3_human",   # latrophilin 3 (ADGRL3)
    "C": "grm5_human",     # mGlu5
    "F": "smo_human",      # Smoothened
    "T2": "t2r46_human",   # TAS2R46
}


def class_of(receptor_class):
    """GPCRdb class name -> template key: 'Class A (Rhodopsin)' -> 'A', 'Class B1 (Secretin)' -> 'B1'."""
    for key in REFERENCE:
        if receptor_class.startswith(f"Class {key} "):
            return key
    raise ValueError(f"no template for {receptor_class!r} yet")

SITE_DATA = Path(__file__).parent.parent / "site" / "data"


def plain(name):
    """GPCRdb names contain HTML ('&beta;<sub>2</sub>'); keep text only."""
    return html.unescape(re.sub(r"<[^>]+>", "", name))


def build(entry_name):
    protein = sources.gpcrdb_protein(entry_name)
    residues = sources.gpcrdb_residues(entry_name)
    model_path, af = sources.alphafold_model(protein["accession"])

    # pLDDT sits in the B-factor column of AlphaFold models.
    chain = gemmi.read_structure(str(model_path))[0][0]
    plddt = {r.seqid.num: round(r["CA"][0].b_iso, 1) for r in chain}

    # Move the model onto the class reference so receptors can be overlaid in 3D.
    gpcr_class = class_of(protein["receptor_class"])
    model_name = f"{af['modelEntityId']}.cif"
    aligned = SITE_DATA / "models" / model_name
    aligned.parent.mkdir(parents=True, exist_ok=True)
    ref_entry = REFERENCE[gpcr_class]
    ref_model, _ = sources.alphafold_model(sources.gpcrdb_protein(ref_entry)["accession"])
    rmsd, n = superpose(model_path, residues, ca_coords(ref_model),
                        sources.gpcrdb_residues(ref_entry), aligned)
    print(f"superposed on {ref_entry}: RMSD {rmsd:.2f} A over {n} TM C-alphas")

    # Camera for the 3D viewer: membrane horizontal, outside up, TM1 on the
    # left and TM7 on the right, like the snake plot.
    _, _, x50, normal, centre, bundle = membrane_frame(residues, ca_coords(aligned))
    right = x50[7] - x50[1]
    right -= np.dot(right, normal) * normal
    right /= np.linalg.norm(right)
    view = {"target": np.round(bundle, 2).tolist(), "up": np.round(normal, 4).tolist(),
            "toward_camera": np.round(np.cross(right, normal), 4).tolist()}

    out_residues = []
    for r in residues:
        g = parse(r["display_generic_number"])
        out_residues.append({
            "seq": r["sequence_number"],
            "aa": r["amino_acid"],
            "segment": r["protein_segment"],
            "generic": g.key if g else None,           # template key, e.g. "3x50"
            "display": r["display_generic_number"],    # e.g. "3.50x50"
            "plddt": plddt.get(r["sequence_number"]),
        })

    data = {
        "entry_name": entry_name,
        "name": plain(protein["name"]),
        "gene": (protein.get("genes") or [""])[0],
        "accession": protein["accession"],
        "class": protein["receptor_class"],
        "class_template": gpcr_class,
        "family": protein["family"],
        "model": {"file": f"data/models/{model_name}", "id": af["modelEntityId"],
                  "version": af["latestVersion"], "mean_plddt": af["globalMetricValue"], "view": view,
                  "superposed_on": ref_entry, "rmsd_A": round(rmsd, 2)},
        "residues": out_residues,
    }
    (SITE_DATA / "receptors").mkdir(parents=True, exist_ok=True)
    out = SITE_DATA / "receptors" / f"{entry_name}.json"
    out.write_text(json.dumps(data, indent=1))
    print(f"wrote {out} ({len(out_residues)} residues) and {model_name}")
    return data


if __name__ == "__main__":
    build(sys.argv[1] if len(sys.argv) > 1 else "adrb2_human")
