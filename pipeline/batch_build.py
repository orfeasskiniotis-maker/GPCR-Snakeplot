"""Build every human receptor of a class, plus the index the home page searches.

Writes:
  site/data/receptors/<entry>.json and site/data/models/<AF id>.cif  (one per receptor)
  site/data/index.json   [{entry, name, gene, accession, family, ligand_type, mean_plddt}, ...]
  pipeline/cache/failed_<class>.json   receptors that could not be built, and why

Usage:  python batch_build.py A      (or B1, B2, C, F, T2)
Takes a while the first time (about 290 receptors, ~120 MB of AlphaFold models);
everything is cached, so re-runs are quick.
"""
import json
import sys
import time
import traceback

import sources
import build_receptor

# GPCRdb family slug of each class. Olfactory (O1, O2), vomeronasal and
# unclassified receptors are left out, as in the roadmap.
CLASS_SLUG = {"A": "001", "B1": "002", "B2": "003", "C": "004", "F": "006", "T2": "009"}
SPECIES = "Homo%20sapiens"


def receptors_in(slug):
    """GPCRdb list of human receptors under a family slug."""
    return json.loads(sources._get(f"{sources.GPCRDB}/proteinfamily/proteins/{slug}/{SPECIES}/",
                                   f"family_{slug}_human.json"))


def ligand_types(slug):
    """{entry_name: 'Aminergic receptors', ...} from the first level below the class."""
    children = json.loads(sources._get(f"{sources.GPCRDB}/proteinfamily/children/{slug}/",
                                       f"family_{slug}_children.json"))
    out = {}
    for child in children:
        for p in receptors_in(child["slug"]):
            out[p["entry_name"]] = child["name"]
    return out


def main(gpcr_class):
    slug = CLASS_SLUG[gpcr_class]
    proteins = receptors_in(slug)
    types = ligand_types(slug)
    # The list already holds what the per-receptor protein call returns: cache it.
    for p in proteins:
        path = sources.CACHE / f"{p['entry_name']}_protein.json"
        if not path.exists():
            path.write_text(json.dumps(p))

    index, failed = [], {}
    start = time.time()
    for i, p in enumerate(proteins, 1):
        entry = p["entry_name"]
        try:
            data = build_receptor.build(entry)
            index.append({
                "entry": entry, "name": data["name"], "gene": data["gene"],
                "accession": data["accession"], "family": data["family"],
                "ligand_type": types.get(entry, ""), "class": gpcr_class,
                "mean_plddt": data["model"]["mean_plddt"],
            })
        except Exception as e:  # keep going; report at the end
            failed[entry] = f"{type(e).__name__}: {e}"
            traceback.print_exc(limit=1)
        print(f"[{i}/{len(proteins)}] {entry} ({time.time() - start:.0f}s)", flush=True)

    # Keep the other classes already in the index; replace this class's entries.
    out = build_receptor.SITE_DATA / "index.json"
    if out.exists():
        index += [r for r in json.loads(out.read_text()) if r.get("class") != gpcr_class]
    order = list(CLASS_SLUG)
    index.sort(key=lambda r: (order.index(r["class"]), r["ligand_type"], r["family"], r["gene"]))
    out.write_text(json.dumps(index, indent=1, ensure_ascii=False))
    (sources.CACHE / f"failed_{gpcr_class}.json").write_text(json.dumps(failed, indent=1))
    print(f"\nbuilt {len(index)} receptors, {len(failed)} failed -> {out}")
    for entry, why in failed.items():
        print(f"  {entry}: {why}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "A")
