#!/usr/bin/env python3
"""Register the three representative MERFISH slice pairs with DET."""

import argparse

import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from sklearn.decomposition import PCA

from _common import ROOT, run


PAIRS = {
    "anterior": ("008", "009", [[-.189903504292, -.981802759752, 10.3056128428],
                                  [.981802759752, -.189903504292, 2.80850849237]]),
    "middle":   ("087", "088", [[-.752838094095, -.658205746009, 3.44893227927],
                                  [.658205746009, -.752838094095, -.795489731486]]),
    "posterior": ("140", "141", [[-.813221626944, -.581954109420, 9.6339991279],
                                   [.581954109420, -.813221626944, 10.9819331524]]),
}
DATA = ROOT / "data/merfish"


def load(number):
    name = f"Zhuang-ABCA-1.{number}"
    return (np.load(DATA / f"{name}-coords.npy"),
            sparse.load_npz(DATA / f"{name}-expression.npz"))


def register(region):
    source_id, target_id, transform = PAIRS[region]
    source, fy = load(source_id)
    target, fx = load(target_id)
    values = PCA(50, svd_solver="randomized", random_state=0).fit_transform(
        sparse.vstack((fy, fx)).toarray())
    fy, fx = np.split(values, [len(source)])
    transform = np.asarray(transform)
    source = source @ transform[:, :2].T + transform[:, 2]

    similarity = run(target, source, fx, fy,
                     "-Tsr -w.01 -g30 -j5 -p -e.5 -DB,2000,.03,.1 "
                     "-c1e-6 -n500 -r1")
    result = run(target, similarity, fx, fy,
                 "-Tsrn -ux -w.1 -b2 -l2 -g.1 -j2 -J300 -K100 -p -e.5 "
                 "-DB,5000,.03,.1 -c1e-6 -n500 -r1")
    return target, source, similarity, result


def show(results):
    fig, axes = plt.subplots(len(results), 3, figsize=(11, 3.6 * len(results)),
                             squeeze=False)
    for row, (region, (target, *sources)) in enumerate(results.items()):
        for ax, source, title in zip(axes[row], sources,
                                     ("Input", "Similarity DET", "Nonrigid DET")):
            ax.scatter(*target.T, s=.15, label="Target", rasterized=True)
            ax.scatter(*source.T, s=.15, color="tab:red", label="Source",
                       rasterized=True)
            ax.set(aspect="equal")
            if row == 0:
                ax.set_title(title)
            ax.axis("off")
        axes[row, 0].text(-.04, .5, region.title(), rotation=90,
                          transform=axes[row, 0].transAxes, va="center")
    axes[0, 0].legend(markerscale=12)
    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("region", nargs="?", default="all",
                        choices=("all", *PAIRS))
    choice = parser.parse_args().region
    regions = PAIRS if choice == "all" else (choice,)
    show({region: register(region) for region in regions})
