# MOSTA data

`mosta.npz` contains the E14.5 source and E15.5 target used by the accepted
DET paper: spatial coordinates, anatomical-domain annotations, and 50 shared
principal components fitted jointly to the two slices.

The file is distributed separately because it is generated from the large
MOSTA `.h5ad` files and is not suitable for Git storage. Extract the MOSTA
demo-data archive here before running `demo-python/det/mosta.py`.
