#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Adsorption‐Site Generator for Prebuilt HEA Slabs

1. Aligns each slab so its surface normal points along +z.
2. Screens the top‐percentile of z‐coordinates to pick candidate surface atoms.
3. Projects candidates into the xy‐plane and extracts the 2D convex hull vertices.
4. Generates top, bridge, and hollow adsorption sites above those hull atoms.
5. Writes one .vasp file per site into the output directory.

All parameters (input/output dirs, screening percentile, distances, etc.) are
configurable via command line.
"""

import os
import argparse
import logging

import numpy as np
from scipy.spatial import ConvexHull
from ase.io import read, write
from ase import Atom
from ase.geometry import find_mic

# Configure logging
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

def align_slab_to_z(slab):
    """Rotate the slab so that its surface normal (a × b) aligns with +z."""
    cell = slab.get_cell()
    n = np.cross(cell[0], cell[1])
    n /= np.linalg.norm(n)
    if n[2] < 0:
        n = -n
    target = np.array([0,0,1])
    cosang = np.clip(np.dot(n, target), -1.0, 1.0)
    angle = np.arccos(cosang)
    if angle > 1e-6:
        axis = np.cross(n, target)
        axis /= np.linalg.norm(axis)
        slab.rotate(angle, axis, center='COM', rotate_cell=True)
        logging.info(f"Rotated slab by {np.degrees(angle):.2f}°")
    return slab

def identify_surface_atoms_via_hull(slab, percentile=70):
    """
    1) Compute z_threshold = percentile(z_i).
    2) Select atom indices with z_i >= z_threshold.
    3) Project their (x,y) coords and compute 2D convex hull.
    4) Return those atom indices on the hull.
    """
    pos = slab.get_positions()
    zs = pos[:,2]
    z_th = np.percentile(zs, percentile)
    cand = np.where(zs >= z_th)[0]
    logging.info(f"{percentile}th percentile z = {z_th:.3f} Å → {len(cand)} candidates")
    if len(cand) < 3:
        logging.warning("Too few candidates; using all atoms as surface.")
        return cand.tolist()
    xy = pos[cand][:,:2]
    hull = ConvexHull(xy)
    surf = cand[hull.vertices]
    logging.info(f"{len(surf)} surface atoms found via convex hull")
    return surf.tolist()

def generate_adsorption_sites(slab, surface_atoms,
                              ads_dist=1.8, dist_thr=3.0, margin=0.2):
    """
    Generate:
      - Top sites: each surface atom + ẑ * ads_dist
      - Bridge: midpoint of any two surface atoms within dist_thr + ẑ * ads_dist
      - Hollow: centroid of any three surface atoms pairwise within dist_thr + ẑ * ads_dist
    Keep only those with z > base_z + margin.
    """
    pos = slab.get_positions()
    cell, pbc = slab.get_cell(), slab.get_pbc()
    normal = np.array([0.0, 0.0, 1.0])
    sites = []

    # Top
    for idx in surface_atoms:
        bz = pos[idx,2]
        site = pos[idx] + normal*ads_dist
        if site[2] > bz + margin:
            sites.append((f"Top_{idx}", site))

    # Bridge
    for i, idx1 in enumerate(surface_atoms):
        for idx2 in surface_atoms[i+1:]:
            vec, d = find_mic(pos[idx2]-pos[idx1], cell, pbc)
            if d <= dist_thr:
                midpoint = pos[idx1] + vec/2
                site = midpoint + normal*ads_dist
                if site[2] > max(pos[idx1,2],pos[idx2,2]) + margin:
                    sites.append((f"Bridge_{idx1}_{idx2}", site))

    # Hollow
    n = len(surface_atoms)
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                idx1, idx2, idx3 = surface_atoms[i], surface_atoms[j], surface_atoms[k]
                v12, d12 = find_mic(pos[idx2]-pos[idx1], cell, pbc)
                v13, d13 = find_mic(pos[idx3]-pos[idx1], cell, pbc)
                v23, d23 = find_mic(pos[idx3]-pos[idx2], cell, pbc)
                if d12<=dist_thr and d13<=dist_thr and d23<=dist_thr:
                    centroid = pos[idx1] + (v12 + v13)/3
                    site = centroid + normal*ads_dist
                    if site[2] > max(pos[idx1,2],pos[idx2,2],pos[idx3,2]) + margin:
                        sites.append((f"Hollow_{idx1}_{idx2}_{idx3}", site))

    logging.info(f"{len(sites)} adsorption sites generated")
    return sites

def add_adsorbate_and_write(slab, label, position, base, out_dir, ads):
    """Append Atom(ads) at position, then write a VASP POSCAR file."""
    slab2 = slab.copy()
    slab2 += Atom(ads, position=position)
    fname = f"{base}_{label}.vasp"
    write(os.path.join(out_dir, fname), slab2, format='vasp')
    logging.info(f"  → Wrote {fname}")

def main():
    p = argparse.ArgumentParser(
        description="Add adsorption sites to existing HEA slabs via convex hull screening"
    )
    p.add_argument("--input_dir",  required=True,
                   help="Directory containing slab (.vasp/.xyz) files")
    p.add_argument("--output_dir", required=True,
                   help="Directory to write slab+adsorbate files")
    p.add_argument("--percentile", type=float, default=70,
                   help="z‐coordinate percentile for candidate screening")
    p.add_argument("--ads_dist",   type=float, default=1.8,
                   help="Distance above surface for adsorbate (Å)")
    p.add_argument("--dist_thr",   type=float, default=3.0,
                   help="Max interatomic distance for bridge/hollow (Å)")
    p.add_argument("--margin",     type=float, default=0.2,
                   help="Height margin to ensure above the surface (Å)")
    p.add_argument("--adsorbate",  type=str,   default="H",
                   help="Element symbol of the adsorbate")
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    for fname in os.listdir(args.input_dir):
        if not (fname.endswith(".vasp") or fname.endswith(".xyz")):
            continue
        path = os.path.join(args.input_dir, fname)
        try:
            slab = read(path)
        except Exception as e:
            logging.error(f"Failed to read {fname}: {e}")
            continue

        slab = align_slab_to_z(slab)
        surf = identify_surface_atoms_via_hull(slab, percentile=args.percentile)
        if not surf:
            logging.warning(f"No surface atoms for {fname}")
            continue

        sites = generate_adsorption_sites(
            slab, surf,
            ads_dist=args.ads_dist,
            dist_thr=args.dist_thr,
            margin=args.margin
        )

        base = os.path.splitext(fname)[0]
        for label, pos in sites:
            add_adsorbate_and_write(
                slab, label, pos, base,
                args.output_dir, args.adsorbate
            )

if __name__ == "__main__":
    main()
