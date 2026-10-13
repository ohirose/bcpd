# Python demos

Build BCPD, install the Python packages, and run a named preset:

```sh
make
python -m pip install -r demo-python/requirements.txt
python demo-python/bcpd/nonrigid.py bunny-a
python demo-python/bcpd/plusplus.py dragon
python demo-python/bcpd/rigid.py chef-a
python demo-python/bcpd/reconstruction.py
python demo-python/bcpd/hierarchy.py female
python demo-python/bcpd/hierarchy.py male
python demo-python/bcpd/transfer.py a
python demo-python/bcpd/gbcpd.py face-1
python demo-python/bcpd/gbcpd_plusplus.py armadillo
python demo-python/dld/face.py 1 1
python demo-python/dld/hand.py 1
```

DET demos use the same executable and temporary-file workflow:

```sh
python demo-python/det/audio.py
python demo-python/det/image.py gray
python demo-python/det/image.py color
python demo-python/det/shape.py localization-1
python demo-python/det/shape.py localization-2
python demo-python/det/shape.py registration
python demo-python/det/merfish.py
python demo-python/det/mosta.py
```

The MERFISH demo runs the paper figure's representative anterior, middle, and
posterior slice pairs. Pass `anterior`, `middle`, or `posterior` to run only one
pair. Unlike the paper experiment, the compact demo fits PCA independently to
each pair.

The MOSTA demo reproduces the paper's DET-only E14.5-to-E15.5 registration
with the same shared 50-component PCA representation and two registration
stages. Its separately distributed data archive must first be extracted into
`data/mosta`.

Visualization is selected automatically: Matplotlib displays 2D point sets,
while Open3D displays full-resolution 3D point sets and meshes. Mesh targets
are blue points and source/deformed shapes are black wireframes, avoiding
depth-buffer artifacts when the two surfaces coincide. Static 3D results open
the before and after views simultaneously; closing either window closes both.
Every trajectory demo opens this comparison after its animation window closes.

The nonrigid presets are `armadillo-a/b`, `bunny-a/b`, `dragon-a/b`,
`face-a/b`, `fish-a/b`, and `monkey-a/b`. The BCPD++ presets are
`armadillo`, `asian-dragon`, `dragon`, and `lucy`. Rigid presets combine
`apartment`, `chef`, `para`, `stairs`, or `trex` with `a` through `i`.
GBCPD presets combine `face`, `female`, or `male` with `1` or `2`; GBCPD++
also provides `armadillo`.

Shape-transfer presets are `a` through `f`, matching the MATLAB shell demos.
The two stages are saved under `demo-python/output/shape-transfer/PRESET`; add
`--no-view` to save them without opening the shaded source, target, and
deformed-source mesh viewers.

DLD face and hand models are trained in memory with leave-one-out validation.
Add `--refine` to fine-tune the DLD result with BCPD.

Generated BCPD files live in a temporary directory and are removed after the
viewer closes. Missing distributed data can be supplied through `data/local`.

## View saved results

Run from the repository root with the Fish files installed under `data/`.
Use the executable directly to keep the output files:

```sh
./bcpd -x data/fish-x.txt -y data/fish-y.txt -w0 -b2 -l2 -g3 -c1e-6 -n500 -r1 -o fish_
python demo-python/viewers/view_results.py data/fish-x.txt data/fish-y.txt fish_y.txt
python demo-python/viewers/view_results.py data/fish-x.txt data/fish-y.txt fish_y.txt --save fish.png
```

The viewer arguments are original target, original source, and final transformed
source. It uses Matplotlib for 2D and sequential Open3D windows for 3D: close
Before to see After. `--save` writes a static 2D/3D image without opening a window;
`--matplotlib` offers an interactive 3D alternative.
On Windows, replace `./bcpd` with `.\win\bcpd.exe`.

## Check a synthetic Fish registration

This 91-point 2D, in-sample teaching example uses known identities to measure
before/after distances. It is not anatomical or held-out validation. The output
folder `fish-control` must not already exist.

```sh
python demo-python/bcpd/landmark_check.py prepare data/fish-x.txt fish-control
./bcpd -x fish-control/target.txt -y fish-control/source.txt -w0 -b2 -l2 -g3 -c1e-6 -n500 -r1 -o fish-control/fit_
python demo-python/bcpd/landmark_check.py measure fish-control fish-control/fit_y.txt
```

Use this control's final `fit_y.txt` without reordering its rows. The identity map
is for evaluation only; it is not supplied to BCPD.
