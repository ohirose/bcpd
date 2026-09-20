#!/usr/bin/env python3
"""Register and compare two audio signals with DET."""

import matplotlib.pyplot as plt
import numpy as np

from _common import data_file, run

RATE = 44100


if __name__ == "__main__":
    target = np.loadtxt(data_file("Ensoniq-ZR-76-Ac-Bass-2-C2.txt"))
    source = np.loadtxt(data_file("Alesis-Sanctuary-QCard-AcoustcBas-C2.txt"))
    x = np.arange(1, len(target) + 1)[:, None]
    y = np.arange(1, len(source) + 1)[:, None]
    result = run(x, y, target, source,
                 "-DB,2000,.05 -A -e.4 -l10 -b1 -w.1 -g1 -r1 "
                 "-ux -Ux -n1000 -c1e-6")

    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharex=True, sharey=True)
    for ax, time, title in zip(axes, (y[:, 0], result),
                               ("Before Registration", "After Registration")):
        ax.plot(time / RATE, source, "r", lw=.5, label="Source")
        ax.plot(x[:, 0] / RATE, target, "b", lw=.5, label="Target")
        ax.set(xlim=(0, 2.5), ylim=(-1, 1), title=title,
               xlabel="Time (s)", ylabel="Amplitude")
        ax.grid()
    axes[0].legend()
    fig.tight_layout()
    plt.show()
