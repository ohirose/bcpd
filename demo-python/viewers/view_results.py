#!/usr/bin/env python3
"""View final BCPD/DET coordinates; --save needs no graphical display."""

import argparse
from pathlib import Path

import numpy as np


def plot(target, source, result, save=None):
    if save:
        import matplotlib
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    dimension = target.shape[1]
    low = np.minimum.reduce([p.min(0) for p in (target, source, result)])
    high = np.maximum.reduce([p.max(0) for p in (target, source, result)])
    width = max(float((high - low).max()), 1e-9)
    fig = plt.figure(figsize=(12, 5.5))
    for i, (points, title) in enumerate(((source, "Before"), (result, "After"))):
        ax = fig.add_subplot(1, 2, i + 1,
                             projection="3d" if dimension == 3 else None)
        ax.scatter(*target.T, s=2, color="#0072b2", label="Target")
        ax.scatter(*points.T, s=2, color="#d94e36", label="Source")
        ax.set_title(title)
        setters = (ax.set_xlim, ax.set_ylim) if dimension == 2 else (
            ax.set_xlim, ax.set_ylim, ax.set_zlim)
        for setter, center in zip(setters, (low + high) / 2):
            setter(center - .55 * width, center + .55 * width)
        if dimension == 2:
            ax.set_aspect("equal", adjustable="box")
        else:
            ax.set_box_aspect((1, 1, 1))
        ax.legend(markerscale=3)
    fig.tight_layout()
    if save:
        fig.savefig(save, dpi=200)
        print(f"Saved {Path(save).resolve()}")
        plt.close(fig)
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target", help="original X coordinates")
    parser.add_argument("source", help="original Y coordinates")
    parser.add_argument("result", help="final PREFIX_y.txt coordinates")
    parser.add_argument("--save", metavar="PNG", help="save a static comparison")
    parser.add_argument("--matplotlib", action="store_true",
                        help="use Matplotlib rather than Open3D for 3D")
    args = parser.parse_args()
    target, source, result = (np.loadtxt(p, ndmin=2) for p in
                              (args.target, args.source, args.result))
    if (target.shape[1] not in (2, 3) or
            source.shape[1] != target.shape[1] or result.shape != source.shape):
        parser.error("need 2D/3D coordinates and one result row per source row")
    if not all(np.isfinite(p).all() for p in (target, source, result)):
        parser.error("coordinates must be finite")
    if target.shape[1] == 2 or args.save or args.matplotlib:
        plot(target, source, result, args.save)
        return
    try:
        import open3d as o3d
    except ImportError:
        parser.error("Open3D is missing; install requirements or use --matplotlib")
    for points, title in ((source, "Before"), (result, "After")):
        clouds = []
        for vertices, color in ((target, [0, .45, .74]),
                                (points, [.85, .1, .1])):
            cloud = o3d.geometry.PointCloud(o3d.utility.Vector3dVector(vertices))
            cloud.paint_uniform_color(color)
            clouds.append(cloud)
        viewer = o3d.visualization.Visualizer()
        if not viewer.create_window(title, width=900, height=700):
            parser.error("cannot create viewer; try --matplotlib or --save result.png")
        for cloud in clouds:
            viewer.add_geometry(cloud)
        viewer.get_render_option().background_color = np.ones(3)
        viewer.get_render_option().point_size = 2
        viewer.run()
        viewer.destroy_window()


if __name__ == "__main__":
    main()
