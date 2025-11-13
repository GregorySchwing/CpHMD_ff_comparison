#!/usr/bin/env python3
"""
compute_solvent_volume_and_ions.py

Compute solvent-only volume (box minus solute excluded volume) using a voxel grid,
then compute the number of ion pairs required for a target molarity.

Usage:
    python compute_solvent_volume_and_ions.py system.prmtop system.rst7 [target_M] [spacing] [probe]

Defaults:
    target_M = 0.150  (150 mM)
    spacing  = 0.8 Å
    probe    = 1.4 Å

Output:
    - prints diagnostics
    - writes tleap_ion_counts.txt with:
        addions model Na+ <N>
        addions model Cl- <M>
"""

from __future__ import annotations
import sys
import math
import time
from typing import Optional, Tuple, Dict

try:
    import numpy as np
except Exception as e:
    print("ERROR: numpy is required (pip install numpy).", file=sys.stderr)
    raise

try:
    import parmed as pmd
except Exception as e:
    print("ERROR: parmed is required (pip install parmed).", file=sys.stderr)
    raise

try:
    from scipy.spatial import cKDTree as KDTree
    _HAVE_KDTREE = True
except Exception:
    KDTree = None
    _HAVE_KDTREE = False

try:
    from tqdm import tqdm
    _HAVE_TQDM = True
except Exception:
    _HAVE_TQDM = False

AVOGADRO = 6.02214076e23

BONDI_RADII = {
    'H':1.20,'C':1.70,'N':1.55,'O':1.52,'F':1.47,
    'P':1.80,'S':1.80,'CL':1.75,'BR':1.85,'I':1.98,
    'NA':2.27,'K':2.75,'MG':1.73,'ZN':1.39,'CA':2.00
}

def vdw_radius_from_symbol(sym: Optional[str]) -> float:
    if not sym:
        return 1.7
    s = sym.upper()
    if s in BONDI_RADII:
        return BONDI_RADII[s]
    if s[0] in BONDI_RADII:
        return BONDI_RADII[s[0]]
    return 1.7

def get_box_lengths(struct):
    box = getattr(struct, 'box', None)
    if box is None:
        return None
    try:
        if len(box) < 3:
            return None
    except Exception:
        return None
    try:
        a = float(box[0]); b = float(box[1]); c = float(box[2])
    except Exception:
        return None
    if a <= 0 or b <= 0 or c <= 0:
        return None
    return a, b, c

def compute_solute_volume_grid(struct, spacing=0.8, probe=1.4,
                               water_resnames=('WAT','HOH','OPC','TIP3','TIP4','SPCE'),
                               use_kdtree=True):
    box_lengths = get_box_lengths(struct)
    if box_lengths is None:
        raise RuntimeError("Box vectors not found in structure; ensure rst7 contains box.")
    a, b, c = box_lengths
    box_vol = a * b * c

    water_set = {w.upper() for w in water_resnames}

    solute_atoms = []
    water_count = 0
    for res in struct.residues:
        if res.name.upper() in water_set:
            water_count += 1
        else:
            for atom in res.atoms:
                try:
                    x = float(atom.xx)
                    y = float(atom.xy)
                    z = float(atom.xz)
                except Exception:
                    x = float(getattr(atom, 'x', 0.0))
                    y = float(getattr(atom, 'y', 0.0))
                    z = float(getattr(atom, 'z', 0.0))
                try:
                    elem = atom.element.symbol if atom.element else None
                except Exception:
                    elem = None
                solute_atoms.append((x, y, z, elem))

    if not solute_atoms:
        raise RuntimeError("No solute atoms found.")

    coords = np.array([[x,y,z] for (x,y,z,_) in solute_atoms], float)
    radii = np.array([vdw_radius_from_symbol(el) + probe for (_,_,_,el) in solute_atoms], float)

    nx = max(1, int(math.ceil(a / spacing)))
    ny = max(1, int(math.ceil(b / spacing)))
    nz = max(1, int(math.ceil(c / spacing)))

    xs = (np.arange(nx) + 0.5) * spacing
    ys = (np.arange(ny) + 0.5) * spacing
    zs = (np.arange(nz) + 0.5) * spacing

    voxel_vol = spacing**3
    occupied = 0

    have_tree = False
    if _HAVE_KDTREE and use_kdtree:
        try:
            tree = KDTree(coords)
            max_r = float(radii.max())
            have_tree = True
        except Exception:
            have_tree = False

    total_cols = nx * ny
    iterator = ((ix,iy) for ix in range(nx) for iy in range(ny))
    if _HAVE_TQDM:
        iterator = tqdm(((ix,iy) for ix in range(nx) for iy in range(ny)),
                        total=total_cols, desc="Columns")

    for ix, iy in iterator:
        x = xs[ix]
        y = ys[iy]
        pts = np.vstack([np.full_like(zs, x),
                         np.full_like(zs, y),
                         zs]).T
        if have_tree:
            idx_lists = tree.query_ball_point(pts, r=max_r)
            for idxs, p in zip(idx_lists, pts):
                if not idxs:
                    continue
                cand_c = coords[idxs]
                cand_r = radii[idxs]
                d2 = ((cand_c - p)**2).sum(axis=1)
                if (d2 <= (cand_r**2)).any():
                    occupied += 1
        else:
            d2 = np.sum((coords[None,:,:] - pts[:,None,:])**2, axis=2)
            mask = (d2 <= (radii**2)[None,:]).any(axis=1)
            occupied += int(mask.sum())

    solute_vol = occupied * voxel_vol
    return {
        'solute_volume_A3': solute_vol,
        'box_volume_A3': box_vol,
        'water_count': water_count,
        'voxel_count': occupied,
        'voxel_vol': voxel_vol,
        'nx': nx, 'ny': ny, 'nz': nz
    }

def compute_ions_from_solvent_volume(solvent_vol_A3, solute_charge_k, target_M):
    V_L = solvent_vol_A3 * 1e-27
    pairs_needed = target_M * AVOGADRO * V_L

    if solute_charge_k < 0:
        Npos_neutral = int(math.ceil(abs(solute_charge_k)))
    else:
        Npos_neutral = 0
    if solute_charge_k > 0:
        Nneg_neutral = int(math.ceil(solute_charge_k))
    else:
        Nneg_neutral = 0

    neutral_pairs = max(Npos_neutral, Nneg_neutral) / 2
    additional_pairs = pairs_needed - neutral_pairs
    if additional_pairs < 0:
        additional_pairs = 0

    add_pairs_int = int(math.ceil(additional_pairs))

    Npos_total = Npos_neutral + add_pairs_int
    Nneg_total = Nneg_neutral + add_pairs_int

    return {
        'pairs_needed_raw': pairs_needed,
        'pairs_needed_ceil': int(math.ceil(pairs_needed)),
        'Npos_neutral': Npos_neutral,
        'Nneg_neutral': Nneg_neutral,
        'neutral_pairs': neutral_pairs,
        'additional_pairs': additional_pairs,
        'additional_pairs_int': add_pairs_int,
        'Npos_total': Npos_total,
        'Nneg_total': Nneg_total,
        'V_L': V_L
    }

def main(argv):
    if len(argv) < 3:
        print("Usage: python compute_solvent_volume_and_ions.py system.prmtop system.rst7 [target_M] [spacing] [probe]")
        sys.exit(1)

    prmtop = argv[1]
    rst7 = argv[2]
    target_M = float(argv[3]) if len(argv) > 3 else 0.150
    spacing = float(argv[4]) if len(argv) > 4 else 0.8
    probe = float(argv[5]) if len(argv) > 5 else 1.4

    print("Loading structure (ParmEd)...")
    struct = pmd.load_file(prmtop, rst7)

    box = get_box_lengths(struct)
    if box is None:
        print("ERROR: Box vectors missing.")
        sys.exit(2)
    a,b,c = box
    box_vol = a*b*c

    print(f"Box: a={a:.3f} b={b:.3f} c={c:.3f}  | Volume = {box_vol:.3f} Å³")

    water_resnames = ('WAT','HOH','OPC','TIP3','TIP4','SPCE')
    water_count = sum(1 for r in struct.residues if r.name.upper() in {w.upper() for w in water_resnames})
    print(f"Water count: {water_count}")

    k = getattr(struct, 'charge', None)
    if k is None:
        k = sum(a.charge for a in struct.atoms)
    k = float(k)
    print(f"Total solute charge k = {k}")

    print(f"Computing solute excluded volume (spacing={spacing}, probe={probe})...")
    t0 = time.time()
    grid = compute_solute_volume_grid(struct, spacing=spacing, probe=probe)
    t1 = time.time()

    solute_vol = grid['solute_volume_A3']
    solvent_vol = box_vol - solute_vol

    print(f"Voxel runtime: {t1 - t0:.1f} s")
    print(f"Solute volume = {solute_vol:.3f} Å³")
    print(f"Solvent-only volume = {solvent_vol:.3f} Å³")

    ion_info = compute_ions_from_solvent_volume(solvent_vol, k, target_M)

    # ---- ACTUAL MOLARITY HERE ----
    actual_pairs = ion_info['additional_pairs_int'] + ion_info['neutral_pairs']
    actual_M = actual_pairs / (AVOGADRO * ion_info['V_L'])

    print("\nIon calculation:")
    print(f"  Target molarity: {target_M} M")
    print(f"  Actual molarity (from ion counts): {actual_M:.6f} M")
    print(f"  TOTAL Na+: {ion_info['Npos_total']}")
    print(f"  TOTAL Cl-: {ion_info['Nneg_total']}")

    outfn = "tleap_ion_counts.txt"
    with open(outfn, "w") as fh:
        fh.write(f"# Target molarity: {target_M} M\n")
        fh.write(f"# Actual molarity: {actual_M:.6f} M\n")
        fh.write(f"# Solvent-only volume (Å^3): {solvent_vol:.6f}\n")
        fh.write(f"# Box volume (Å^3): {box_vol:.6f}\n")
        fh.write(f"addions model Na+ {ion_info['Npos_total']}\n")
        fh.write(f"addions model Cl- {ion_info['Nneg_total']}\n")

    print(f"\nWrote: {outfn}")
    print("Done.")

if __name__ == "__main__":
    main(sys.argv)

