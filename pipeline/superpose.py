"""Superpose a receptor's AlphaFold model onto the class reference model.

AlphaFold models each come in their own coordinate frame. To show two
receptors on top of each other in 3D, every model is moved onto the
reference receptor of its class (the one the template was built from),
using the C-alpha atoms of TM residues that share a generic number.
This is the same idea as the 2D template, applied in 3D.

The coordinates are rewritten inside the original mmCIF file, so all the
other AlphaFold data in it (pLDDT and so on) is kept.
"""
import gemmi
import numpy as np

from generic_numbers import parse


def tm_generic_to_seq(residues):
    """{'3x50': 131, ...} for TM residues only."""
    out = {}
    for r in residues:
        g = parse(r["display_generic_number"])
        if g and r["protein_segment"].startswith("TM"):
            out[g.key] = r["sequence_number"]
    return out


def kabsch(moving, fixed):
    """Rotation R and translation t that best map moving onto fixed (least squares)."""
    mc, fc = moving.mean(axis=0), fixed.mean(axis=0)
    h = (moving - mc).T @ (fixed - fc)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))   # avoid a mirror image
    r = vt.T @ np.diag([1, 1, d]) @ u.T
    return r, fc - r @ mc


def superpose(model_path, residues, ref_ca, ref_residues, out_path):
    """Write model_path, moved onto the reference, to out_path. Returns (RMSD in Å, n atoms)."""
    from build_template import ca_coords
    ca = ca_coords(model_path)
    a, b = tm_generic_to_seq(residues), tm_generic_to_seq(ref_residues)
    shared = sorted(set(a) & set(b))
    moving = np.array([ca[a[k]] for k in shared])
    fixed = np.array([ref_ca[b[k]] for k in shared])
    r, t = kabsch(moving, fixed)
    rmsd = float(np.sqrt(((moving @ r.T + t - fixed) ** 2).sum(axis=1).mean()))

    doc = gemmi.cif.read(str(model_path))
    block = doc.sole_block()
    table = block.find("_atom_site.", ["Cartn_x", "Cartn_y", "Cartn_z"])
    for row in table:
        xyz = np.array([float(row[0]), float(row[1]), float(row[2])]) @ r.T + t
        row[0], row[1], row[2] = (f"{v:.3f}" for v in xyz)
    doc.write_file(str(out_path))
    return rmsd, len(shared)
