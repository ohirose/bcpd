# Demo
If you are a MATLAB user, numerous examples will be demonstrated
by running the scripts in this directory. If your environment is Mac or Linux,
set the value 'win' to 0 before running the demo scripts. The script
`demoPrepare.sh` automates this procedure.

- Nonrigid/BCPD examples:
  - Fish, Monkey, Bunny, Face, Armadillo, Dragon.
- Nonrigid/BCPD++ examples:
  - Armadillo, Dragon, Asian Dragon, Lucy.
- Rigid/BCPD examples:
  - Chef, Parasaurolophus, T-rex, Apartment, Stairs.
- Shape transfer/BCPD++ examples:
  - Human shapes.

## NOTES
- To run the demo regarding the "shape transfer" with Windows, it is essential to
prepare a UNIX-like environment. Please use MINGW and MSYS systems, for example.
- For Windows, demoPlusPlusLucy.m fails because of memory allocation error,
which originates from the MINGW's 32-bit compilation.
