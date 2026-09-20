#!/usr/bin/env python3
"""Reproduce the DET-only MOSTA E14.5-to-E15.5 comparison."""

import matplotlib.pyplot as plt
import numpy as np

from _common import data_file, run


def main():
    alpha = .5
    data = np.load(data_file("mosta/mosta.npz"))
    source, target = data["source"], data["target"]
    fy, fx = data["source_pca"], data["target_pca"]

    similarity = run(target, source, fx, fy,
                     "-r1 -n200 -c1e-6 -Tsr -w.2 -g1 "
                     "-DB,2000,.01 -p -f.3")
    result = run(target, similarity, fx, fy,
                 "-r1 -ux -n200 -c1e-6 -l1 -b1 -w.1 -g.1 "
                 "-J300 -K150 -p -e.5 -DB,10000,.03,.01,.2")

    labels = dict.fromkeys(data["source_labels"])
    colors = {label: plt.cm.tab20(i / len(labels))
              for i, label in enumerate(labels)}
    source_colors = [colors[x] for x in data["source_labels"]]
    target_colors = [colors.get(x, (.8, .8, .8, 1))
                     for x in data["target_labels"]]
    panels = ((similarity, source_colors, "(a) E14.5 (similarity)"),
              (target, target_colors, "(b) E15.5 (target)"),
              (result, source_colors, "(c) E14.5 (nonrigid)"))

    fig, axes = plt.subplots(1, 5, figsize=(18, 4))
    for ax, (points, color, title) in zip(axes, panels):
        ax.scatter(points[:, 0], -points[:, 1], c=color, s=.1,
                   alpha=alpha, rasterized=True)
        ax.set_title(title)
    for ax, points, title in zip(axes[3:], (similarity, result),
                                 ("(d) Before nonrigid", "(e) After nonrigid")):
        ax.scatter(target[:, 0], -target[:, 1], c=target_colors, s=.1,
                   alpha=alpha, rasterized=True)
        ax.scatter(points[:, 0], -points[:, 1], c=source_colors, s=.1,
                   alpha=alpha, rasterized=True)
        ax.set_title(title)
    for ax in axes:
        ax.set_aspect("equal")
        ax.axis("off")
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
