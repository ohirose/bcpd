
# Bayesian Coherent Point Drift / Bayesian Coherent Point Drift ++

This is an implementation of a non-rigid point matching algorithm, Bayesian coherent point drift (BCPD), with
accelerations based on the Nystrom method and the KD tree search. BCPD combines non-rigid and rigid registration.
Therefore,
(1) BCPD solves non-rigid registration with robustness against target rotation and
(2) BCPD solves rigid registration under an appropriate set of tuning parameters.
The algorithm can further be accelerated using downsampling and interpolation. We call the acceleration scheme BCPD++.
+Algorithm details are available in [Hirose2020a](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8985307)
+and [Hirose2020b](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9290402).

## Table of Contents

1. [Papers](#papers)
2. [Demo](#demo)
3. [Usage](#usage)
    + [Terms and Symbols](#terms-and-symbols)
    + [Input data](#input-data)
4. [Options](#options)
    + [Tuning parameters](#tuning-parameters)
    + [Kernel functions](#kernel-functions)
    + [Acceleration](#acceleration)
    + [Downsampling](#downsampling)
    + [Convergence](#convergence)
    + [Normalization](#normalization)
    + [File output](#file-output)
    + [Terminal output](#terminal-output)
5. [Rigid registration](#rigid-registration)

## Papers

The details of the algorithms are available in the following papers:
- [BCPD++] O. Hirose,
  "[Acceleration of non-rigid point set registration with downsampling and Gaussian process regression](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=9290402)",
  IEEE TPAMI, Dec 2020.
- [BCPD] O. Hirose,
  "[A Bayesian formulation of coherent point drift](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8985307)",
  IEEE TPAMI, Feb 2020.
  - Supplementary Video 1 in the above paper is available [HERE](https://youtu.be/cET6gKAvjw0).
  If the video file cannot be accessed, go to [online-materials](https://ieeexplore.ieee.org/document/8985307/media#media).

## Demo

If you are a MATLAB user, demo codes can be executed in the command window of MATLAB.
- Download the datasets required for demos, which are available
  [HERE](https://www.dropbox.com/s/6kd4uiyt150uyz9/bcpd-demodata20200127.zip?dl=1).
  If you have trouble downloading them, go to [bcpd-dataset](https://github.com/ohirose/bcpd-dataset).
- Decompress and move the datasets into the `data` folder in this software.
- Start MATLAB.
- Go to the `demo` folder in the MATLAB environment.
- Double-click a demo script, e.g., `demoAsianDragon.m`.
- Press the run button in the code editor of MATLAB.

## Usage

Type the following command in the DOS prompt:

` bcpd -x <target: X> -y <source: Y> (+options) `

Brief instructions are printed by typing `./bcpd -v` (or `bcpd -v` for windows) in the terminal window.
The binary file can also be executed using `system` function in MATLAB.
See MATLAB scripts in the `demo` folder regarding the usage of the binary file.

### Terms and symbols

- X: Target point set. The point set corresponding to the reference shape.
- Y: Source point set. The point set to be deformed. The mth point in Y is denoted by ym.
- N: The number of points in the target point set.
- M: The number of points in the source point set.
- D: Dimension of the space in which the source and target point sets are embedded.

### Input data

- 1st argument (specified by `-x`): The target shape represented as a matrix of size N x D.
- 2nd argument (specified by `-y`): The source shape represented as a matrix of size M x D.

Currently, only tab- and comma-separated files are accepted, and the extensions of input files
MUST be `.txt`. If your file is space-delimited, convert it to a tab- or comma-separated file using Excel,
MATLAB, or R, for example. If the file names of target and source point sets are `X.txt` and `Y.txt`,
these arguments can be omitted.

## Options

The tuning parameters and options are listed in the following sections. Default values
will be used if they are not specified.

### Tuning parameters

- `-w [real]`: Omega. Outlier probability in (0,1).
- `-l [real]`: Lambda. Positive. It controls the expected length of displacement vectors.
- `-k [real]`: Kappa. Positive. It controls the randomness of mixing coefficients.
- `-g [real]`: Gamma. Positive. It defines the randomness of the point matching during the early stage of the optimization.

BCPD is a unified framework of non-rigid registration and rigid registration.
If point sets to be registered are smooth surfaces of 3D models, set `-w 0`.
If your target point set is largely rotated, set gamma around
2 to 10, which often contributes to converge a better solution.
If lambda (-l) is sufficiently large, e.g. 1e9, BCPD solves rigid registration problems.
If you would like to solve rigid registration for large point sets, accelerate the algorithm carefully;
see [Rigid registration](#rigid-registration).

### Kernel functions

- `-G [1-5]`: Switch kernel functions. The Gaussian kernel `exp(-||ym-ym'||^2/2*beta^2)` is used unless specified.
  - `-G1` Inverse multiquadric: `(||ym-ym'||^2+beta^2)^(-1/2)`
  - `-G2` Rational quadratic: `1-||ym-ym'||^2/(||ym-ym'||^2+beta^2)`
  - `-G3` Laplace: `exp(-|ym-ym'|/beta)`
  - `-G4` Neural network: see [Williams, Neural computation, 1998] for the definition of the kernel.
- `-b [real(s)]`: The parameter(s) of a kernel function.
  - `-b [real]`: Beta. The parameter of a kernel function except the neural network kernel.
  - `-b [real,real]`: The parameters of the neural network kernel. Do not insert whitespaces before and after comma.

Here, `ym` represents the mth point in Y. Except the neural network kernel, the tuning parameter of
the kernel functions is denoted by beta, which controls the directional correlation of displacement vectors.
If the kernel is Gaussian, the expected length of displacement vectors is controlled by lambda regardless of beta.
Then, the expected length equals to sqrt(D/lambda). For the neural network kernel, the first and second arguments of
the option `-b` specify the standard deviations of the intercept and linear coefficients, respectively.

### Acceleration

- `-A`: Acceleration with default acceleration parameters, i.e., `-K70 -J300 -p -d7 -e0.15 -f0.2`.
- `-K [int]`: #Nystrom samples for computing G.
- `-J [int]`: #Nystrom samples for computing P.
  - `-r [int]`: Random number seed for the Nystrom method. Reproducibility is guaranteed if the same number is specified.
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
are `-d 7`, `-e 0.15` and `-f 0.2`. The computational load will, therefore, be relaxed by specifying
`-d 4` or `-e 0.1`, for example, although the accuracy of the computation decreases.
If J, K, e, and d are not enough, the optimization will become unstable.

### Downsampling

- `-D [char,int,real]`: Changes the number of points. E.g., `-D'B,10000,0.08'`.
  - 1st argument: One of the symbols: [X,Y,B,x,y,b]; x: target; y: source; b: both, upper: voxel grid, lower: inverse density.
  - 2nd argument: The number of points to be extracted by the downsampling.
  - 3rd argument: The parameter of a downsampling technique based on the inverse point distribution.
- `-L [int]`: #Nystrom samples for accelerating interpolation (BCPD++).

The algorithm can also be accelerated by downsampling techniques: i) voxel grid filter with voxel width r,
ii) the inverse point distribution with the radius parameter r, and iii) the random sampling with equivalent
sampling probabilities. The parameter r can be specified as the 3rd argument of `-D`. If r is specified as 0,
sampling scheme iii) is selected. Sampling scheme ii) is slightly accurate but much slower than sampling scheme i).
The algorithmic accelerations and a downsampling technique can be combined. The resulting registered shape with
interpolation is output to the file with suffix `y.interpolated.txt`. The numbers of points to be downsampled for
target and source point sets can be different; specify the `-D` option twice, e.g., `-D'X,6000,0.08' -D'Y,5000,0.05'`.

### Convergence

- `-c [real]`: Convergence tolerance.
- `-n [int ]`: The maximum number of VB loops.
- `-N [int ]`: The minimum number of VB loops.

The default value of the convergence tolerance is `1e-4`. If your point sets are smooth
surfaces with moderate numbers of points, specify `-c 1e-5` or `-c 1e-6`.

### Normalization

- `-u [char]`: Chooses a normalization option by specifying the argument of the option, e.g., `-ux`.
  - `e`: Each of X and Y is normalized separately (default).
  - `x`: X and Y are normalized using the location and the scale of X.
  - `y`: X and Y are normalized using the location and the scale of Y.
  - `n` : Normalization is skipped.

### File output

- `-o [string]`: Prefix of file names to be output.
- `-s [string]`: Save variables by specifying them as the argument of the option, e.g., `-sYP`.
  - `y`: Resulting deformed shape (=y).
  - `x`: Target shape with alignment (=x).
  - `u`: Deformed shape without similarity transformation (=u).
  - `v`: Displacement vector (=v).
  - `c`: non-outlier labels (=c).
  - `e`: matched points (=e).
  - `a`: Mixing coefficients (=alpha).
  - `P`: Nonzero matching probabilities (=P).
  - `T`: Similarity transformation (=s,R,t).
  - `Y`: Optimization trajectory.
  - `t`: Computing time (real/cpu) and sigma for each loop.
  - `A`: All of the above.

The resulting deformed shape y will be output without `-s` option. Shape x is roughly the same
as y if two point sets are successfully registered. If at least one of `u`,`v`, and `T` is
specified as an argument of `-s`, normalized X and Y before optimization, which are used as
inputs of BCPD, will be output besides the variables. If `Y` is specified as an argument of
`-s`, the optimization trajectory will be saved to the binary file `.optpath.bin`.
The trajectory can be viewed using the following MATLAB scripts, `optpath.m` for 2D data and
`optpath3.m` for 3D data. Saving a trajectory is memory-inefficient. Disable it if both N and M
are more than several hundreds of thousands. If `P` is specified as an argument of `-s`,
nonzero elements of matching probability P will be output. If the optimization is not converged,
output of P might become time-consuming.

### Terminal output

- `-q`: Quiet mode. Print nothing.
- `-h`: History mode. Status information for each loop will not be cleared if specified.
- `-v`: Print the version and the simple instruction of this software.

## Rigid registration

BCPD solves rigid registration problems if lambda is sufficiently large, e.g. 1e9. To stabilize
the registration performance of the rigid registration, accelerate the algorithm carefully.
For example, use the following option:

- `-l1e9 -w0.1 -J300 -K70 -p -e0.3 -f0.3 -g3 -DB,2000,0.08 -sY`.

Otherwise, the computation will be unstable. If two point sets are roughly registered,
it is a good choice to use `-g0.1 -ux` instead of `-g3`. Do not output P, i.e., specify neither
`-sP` nor `-sA` because the number of nonzero elements in P will be enormous.

