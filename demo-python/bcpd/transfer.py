#!/usr/bin/env python3
"""Transfer a source mesh to a target shape in two BCPD++ stages."""

import argparse

import numpy as np

from _common import ROOT, register
from meshpath import transfer

DATA = ROOT / "demo-matlab/shapeTransfer"
CASES = {
    "a": ("01", "42", 70, ""), "b": ("01", "20", 70, ""),
    "c": ("09", "20", 70, ""), "d": ("01", "09", 70, ""),
    "e": ("20", "13", 100, ""), "f": ("01", "13", 70, "-L100"),
}


def load_obj(name):
    rows = [line.split() for line in (DATA / f"{name}.obj").read_text().splitlines()]
    vertices = np.array([row[1:4] for row in rows if row[0] == "v"], float)
    faces = np.array([[int(i.split("/")[0]) - 1 for i in row[1:]]
                      for row in rows if row[0] == "f"])
    return vertices, faces


def save_obj(path, vertices, faces):
    np.savetxt(path, vertices, fmt="v %.16g %.16g %.16g")
    with path.open("ab") as file:
        np.savetxt(file, faces + 1, fmt="f" + " %d" * faces.shape[1])


def triangles(faces):
    return np.concatenate([faces[:, (0, i, i + 1)]
                           for i in range(1, faces.shape[1] - 1)])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=CASES)
    parser.add_argument("--no-view", action="store_true")
    args = parser.parse_args()

    target, source, rank, extra = CASES[args.demo]
    x, xfaces = load_obj(target)
    y, yfaces = load_obj(source)
    base = f"-J300 -K{rank} -p -r1 {extra} -c1e-6 -n1000 -l50 -DB,4000,0.08"
    stages = (register(x, y, f"{base} -g10 -b2"),)
    stages += (register(x, stages[0], f"{base} -g.1 -b1.2 -ux"),)

    folder = ROOT / "demo-python/output/shape-transfer" / args.demo
    folder.mkdir(parents=True, exist_ok=True)
    for number, vertices in enumerate(stages, 1):
        output = folder / f"transferV{number}_y.obj"
        save_obj(output, vertices, yfaces)
        print(f"Saved {output}")
    if not args.no_view:
        transfer(x, triangles(xfaces), y, triangles(yfaces), stages[-1])
