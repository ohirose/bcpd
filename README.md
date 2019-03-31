
# Bayesian Coherent Point Drift

This is an implementation of a point matching algorithm, Bayesian coherent point drift (BCPD), with
accelerations based on the Nystrom method and the KD tree search. The BCPD is an extension of the coherent
point drift (CPD) [Myronenko and Song, 2010], and the main difference between them lies in their
formulations. The BCPD unifies non-rigid transformation and similarity transformation, and thereby,
the BCPD is often robust against the rotation of target shape. The details of the algorithm are
available at XXX. Currently, we distribute the windows version only.

## Demo

If you are a MATLAB user, demo codes can be executed in the command window of MATLAB.

- Start MATLAB.
- Go to the `demo` folder in the MATLAB environment.
- Double-click a demo script, e.g., `demoBcpdFishA.m`.
- Press the run button in the code editor of MATLAB.

## Usage

Type the following command in the DOS prompt:

` bcpd -x <target: X> -y <source: Y> (+options) `

Brief instructions are printed by typing `bcpd -v` in the terminal window.
The binary file can also be executed using `system` function in MATLAB.
See MATLAB scripts in the `demo` folder about the usage of the binary file.
The binary file was compiled by GCC included in the 32-bit version of the MinGW system,

### Terms and symbols

- X: Target point set. The point set corresponding to the reference shape.
- Y: Source point set. The point set to be deformed.
- N: The number of points in the target point set.
- M: The number of points in the source point set.
- D: Dimension of the space in which the source and target point sets are embedded.

### Input data (required)

- 1st argument (specified by `-x`): The target shape represented as a matrix of size N x D.
- 2nd argument (specified by `-y`): The source shape represented as a matrix of size M x D.

Currently, only tab- and comma-separated files are accepted, and the extensions of input files
MUST be `.txt`. If your file is space-delimited, convert it to tab- or comma-separated using Excel,
MATLAB or R, for example. If the file names of target and source point sets are `X.txt` and `Y.txt`,
these arguments can be omitted.

## Options

The turning parameters and options are listed in the following sections. Default values
will be used if they are not specified.

### Tuning parameters

- `-w [real]`: Omega. Outlier probability in (0,1).
- `-l [real]`: Lambda. Positive. It controls the expected length of displacement vectors.
- `-k [real]`: Kappa. Positive. It controls the randomness of mixing coefficients.
- `-g [real]`: Gamma. Positive. It defines the randomness of the point matching during the early stage of the optimization.

The expected length that is controlled by lambda equals to`sqrt((1/lambda)*D)`.
The BCPD is a unified framework of non-rigid transformation and rigid transformation.
If lambda (-l) is set to a large value, the BCPD solves rigid registration problems.
If point sets to be registered are smooth surfaces of 3D models, set `-w 0`.
If your target point set is largely rotated, set gamma around
2 to 10, which often contributes to converge a better solution.

### Kernel functions

- `-G [1-2]`: Switch kernel functions. The Gaussian kernel `exp(-||ym-ym'||^2/2*beta^2)` is used unless specified.
  - `-G1` Inverse multiquadric: `sqrt(||ym-ym'||^2+beta^2)`
  - `-G2` Rational quadratic: `1-||ym-ym'||/(||ym-ym'||+beta^2)`
- `-b [real]`: Beta. Tuning parameter of the kernel functions.

The kernel and its tuning parameter, beta, controls the directional correlation of displacement
vectors. If the kernel is Gaussian, the expected length of displacement vectors is controlled
by lambda regardless of beta. The expected length equals to sqrt(D/lambda).

### Acceleration

- `-K [int]`: #Nystrom samples for computing G.
- `-J [int]`: #Nystrom samples for computing P.
- `-p`: KD-tree search is turned on if specified. The following options fine-tune the KD tree search.
  - `-d [real]`: Scale factor of sigma that defines areas to search for neighbors.
  - `-e [real]`: Maximum radius to search for neighbors.
  - `-f [real]`: The value of sigma at which the KD tree search is turned on.

The Nystrom method accelerates the execution by a random sampling scheme.
It usually runs faster than the direct computation does if M and N are moderately large
and the number of points to be sampled is set to much smaller than both N and M.
If N and M are larger than several thousand, specify `-J 300 -K 80 -p`, for example.
Then, the computation will be faster without sacrificing the registration accuracy.
Also, we note that N and M are more than several hundreds of thousands, the optimization might
get slow especially near convergence even if the options `-J`, `K`, and `-p` are activated.
The default settings of the scale factor and the maximum radius regarding the KD tree search
are `-d 7` and `-e 0.15`. The computational load will, therefore, be relaxed by specifying
`-d 4` or `-e 0.1`, for example, although the accuracy of the computation decreases.
If J, K, e, and d are not enough, the optimization will become unstable.

### Convergence

- `-c [real]`: Convergence tolerance.
- `-n [int ]`: The maximum number of VB loops.

The default value of the convergence tolerance is `1e-4`. If your point sets are smooth
surfaces with moderate numbers of points, specify `-c 1e-5` or `-c 1e-6`.

### File output

- `-o [char*]`: Prefix of file names to be output.
- `-s [char*]`: Save variables by specifying them as the argument of the option, e.g., `-sYP`.
  - `u`: Normalized deformed shape (=uhat).
  - `x`: Target shape with alignment (=xhat).
  - `a`: Mixing coefficients (=alpha).
  - `P`: Matching probability (=P).
  - `Y`: Optimization trajectory.
  - `A`: All of the above.

The resulting deformed shape yhat will be output without `-s` option. Shape xhat is roughly
the same as yhat if two point sets are successfully registered. If `Y` is specified as an
argument of `-s`, the optimization trajectory will be saved to the binary file `.optpath.bin`.
The trajectory can be viewed using the following MATLAB scripts, `optpath.m` for 2D data and
`optpath3.m` for 3D data. Saving a trajectory is memory-inefficient. Disable it if both N and M
are more than several hundreds of thousands. If `P` is specified as an argument of `-s`,
nonzero elements of matching probability P will be output. For each line of the output file
regarding P, the 1st, 2nd, and 3rd columns represent the index of a point in Y, the index of
a point in X, and their matching probability, respectively.

### Terminal output

- `-q`: Quiet mode. Print nothing.
- `-v`: Print the version and the simple instruction of this software.
- `-h`: Status information for each loop will not be cleared if specified.

