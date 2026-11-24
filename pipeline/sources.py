"""Download helpers for the public data sources the portal uses.

Every function caches its download in pipeline/cache/ so re-running the
pipeline does not hit the APIs again. Delete the cache folder to refresh.
"""
import json
import os
from pathlib import Path

import requests

CACHE = Path(os.environ.get("GPCR_CACHE", Path(__file__).parent / "cache"))
# AlphaFold DB refuses requests that carry no browser-like User-Agent.
HEADERS = {"User-Agent": "Mozilla/5.0 (gpcr-snakeplot pipeline)"}

GPCRDB = "https://gpcrdb.org/services"
AFDB = "https://alphafold.ebi.ac.uk"


def _get(url, cache_name, binary=False):
    """GET a URL once and keep the result in the cache folder."""
    CACHE.mkdir(exist_ok=True)
    path = CACHE / cache_name
    if not path.exists():
        r = requests.get(url, headers=HEADERS, timeout=120)
        r.raise_for_status()
        path.write_bytes(r.content)
    return path.read_bytes() if binary else path.read_text()


def gpcrdb_protein(entry_name):
    """Receptor metadata: name, accession, class, family, sequence."""
    return json.loads(_get(f"{GPCRDB}/protein/{entry_name}/", f"{entry_name}_protein.json"))


def gpcrdb_residues(entry_name):
    """One record per residue: sequence number, amino acid, segment, generic number."""
    return json.loads(_get(f"{GPCRDB}/residues/extended/{entry_name}/", f"{entry_name}_residues.json"))


def alphafold_entry(accession):
    """AlphaFold DB metadata for a UniProt accession (includes the model file URLs)."""
    return json.loads(_get(f"{AFDB}/api/prediction/{accession}", f"{accession}_afdb.json"))[0]


def alphafold_model(accession):
    """Path to the cached AlphaFold mmCIF model for a UniProt accession."""
    entry = alphafold_entry(accession)
    name = f"{entry['modelEntityId']}.cif"
    _get(entry["cifUrl"], name, binary=True)
    return CACHE / name, entry


def gpcrdb_structures():
    """Every experimental GPCR structure GPCRdb knows: PDB code, receptor, class, state, resolution."""
    return json.loads(_get(f"{GPCRDB}/structure/", "gpcrdb_structures.json"))


OPM_PDB = "https://opm-assets.storage.googleapis.com/pdb"


def opm_structure(pdb_code):
    """Path to the OPM copy of a PDB entry, or None if OPM does not have it.

    OPM (Orientations of Proteins in Membranes) rotates and shifts each
    structure so the membrane normal is the z axis and z = 0 is the middle
    of the bilayer. A header line gives half the hydrophobic thickness.
    """
    (CACHE / "opm").mkdir(parents=True, exist_ok=True)
    path = CACHE / "opm" / f"{pdb_code.lower()}.pdb"
    missing = CACHE / "opm" / f"{pdb_code.lower()}.missing"
    if missing.exists():
        return None
    if not path.exists():
        r = requests.get(f"{OPM_PDB}/{pdb_code.lower()}.pdb", headers=HEADERS, timeout=120)
        if r.status_code in (403, 404):
            missing.touch()
            return None
        r.raise_for_status()
        path.write_bytes(r.content)
    return path
