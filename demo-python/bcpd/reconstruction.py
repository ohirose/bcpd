#!/usr/bin/env python3
"""Reconstruct the chef model from 22 point-cloud views."""

import matplotlib.pyplot as plt
import numpy as np

from _common import data_file, rigid_transform, show_cloud

PARENT = (0, *range(1, 16), 3, 8, 1, 8, 4, 2)
OPTIONS = "-w.2 -g1 -J300 -p -f.3 -DB,100000,.04 -c1e-7 -n500 -r1 -Tr -sT"
MAX_POINTS = 300_000


if __name__ == "__main__":
    rotations, translations = {}, {}
    for target in range(2, 23):
        source = PARENT[target - 1]
        rotations[target], translations[target] = rigid_transform(
            f"chef_view{target:03d}.txt", f"chef_view{source:03d}.txt", OPTIONS)

    clouds = []
    for view in range(1, 23):
        points = np.loadtxt(data_file(f"chef_view{view:03d}.txt"))
        while PARENT[view - 1]:
            points = (points - translations[view]) @ rotations[view]
            view = PARENT[view - 1]
        clouds.append(points)

    points = np.vstack(clouds)
    if len(points) > MAX_POINTS:
        points = points[np.random.default_rng(0).choice(
            len(points), MAX_POINTS, replace=False)]
    height = (points[:, 2] - points[:, 2].min()) / np.ptp(points[:, 2])
    show_cloud(points, plt.cm.turbo(height)[:, :3])
