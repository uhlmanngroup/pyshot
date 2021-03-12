#! /usr/bin/env python
import numpy as np
import trimesh
import argparse
import pyshot
import os
import matplotlib.pyplot as plt

"""
A simple CLI to show how to extract SHOT descriptors in python
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser("example_pyshot", description=__doc__)
    parser.add_argument("mesh_file", type=str,
                        help="File containing the mesh")
    parser.add_argument("--radius", type=float, default=100)
    parser.add_argument("--local_rf_radius", type=float, default=None)
    parser.add_argument("--min_neighbors", type=int, default=4)
    parser.add_argument("--n_bins", type=int, default=20)
    parser.add_argument("--double_volumes_sectors", action='store_true')
    parser.add_argument("--use_interpolation", action='store_true')
    parser.add_argument("--use_normalization", action='store_true')

    args = parser.parse_args()

    mesh = trimesh.load(args.mesh_file)

    v = np.array(mesh.vertices)
    f = np.array(mesh.faces)

    local_rf_radius = args.radius if args.local_rf_radius is None else args.local_rf_radius

    shot_descrs = pyshot.get_descriptors(v,
                                         f,
                                         radius=args.radius,
                                         local_rf_radius=local_rf_radius,
                                         min_neighbors=args.min_neighbors,
                                         n_bins=args.n_bins,
                                         double_volumes_sectors=args.double_volumes_sectors,
                                         use_interpolation=args.use_interpolation,
                                         use_normalization=args.use_normalization,
                                         )

    plt.imshow(shot_descrs.T)
    plt.title(f"SHOT descriptors of {os.path.basename(args.mesh_file)} (transposed)")
    plt.show()
