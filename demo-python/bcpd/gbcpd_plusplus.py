#!/usr/bin/env python3
"""Run one of the original GBCPD++ demos."""

from _common import execute

CASES = {
    "armadillo": ("armadillo-g-y.txt", "armadillo-g-x.txt",
                  "armadillo-g-triangles.txt",
                  "-w0 -b1 -l50 -g.1 -J300 -K200 -p -ux -DB,10000,0.02 "
                  "-c1e-6 -n500 -r1 -ux -Ggeodesic,1,{faces}"),
    "face-1": ("face-g-x.txt", "face-g-y.txt", "face-g-triangles.txt",
               "-w0 -b0.7 -l100 -g1 -J300 -K100 -p -ux -DB,5000,0.02 "
               "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
    "face-2": ("face-g-y.txt", "face-g-x.txt", "face-g-triangles.txt",
               "-w0 -b0.7 -l100 -g1 -J300 -K100 -p -ux -DB,5000,0.02 "
               "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
    "female-1": ("female-x.txt", "female-y.txt", "female-triangles.txt",
                 "-w0 -b1.2 -l100 -g3 -J300 -K300 -p -ux -DB,3000,0.02 "
                 "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
    "female-2": ("female-y.txt", "female-x.txt", "female-triangles.txt",
                 "-w0 -b1.2 -l100 -g3 -J300 -K300 -p -ux -DB,3000,0.02 "
                 "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
    "male-1": ("male-x.txt", "male-y.txt", "male-triangles.txt",
               "-w0 -b0.7 -l100 -g1 -J300 -K300 -p -ux -DB,3000,0.02 "
               "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
    "male-2": ("male-y.txt", "male-x.txt", "male-triangles.txt",
               "-w0 -b1.3 -l100 -g3 -J300 -K300 -p -ux -DB,3000,0.02 "
               "-c1e-6 -n500 -r1 -ux -Ggeodesic,.5,{faces}"),
}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("demo", choices=CASES)
    args = parser.parse_args()
    x, y, faces, options = CASES[args.demo]
    execute(x, y, options, faces=faces)
