// Copyright (c) 2018-2020 Osamu Hirose
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include"bcpd.h"
#include"version.h"

#ifdef USE_OPENMP
static int omp=1;
#else
static int omp=0;
#endif

void printUsage(void){
  printf("\n");
  printf(" BCPD version %s (%s). OpenMP: %s.\n",_VERSION_,_DATE_,omp?"Turned on":"Turned off");
  printf(" Copyright (c) Osamu Hirose                                                 \n\n");
  printf(" This software executes the following registration algorithms:              \n"  );
  printf("  o------------------------------------------------------------------------o\n"  );
  printf("  | BCPD    | Bayesian coherent point drift.                               |\n"  );
  printf("  | BCPD++  | A faster BCPD with downsampling and interpolation.           |\n"  );
  printf("  | GBCPD   | Geodesic-based BCPD.                                         |\n"  );
  printf("  | GBCPD++ | A faster GBCPD with downsampling and interpolation.          |\n"  );
  printf("  o------------------------------------------------------------------------o\n\n");
  printf(" USAGE:                                                                     \n"  );
  printf("  o------------------------------------------------------------------------o\n"  );
  printf("  | ./bcpd -x <target point set> -y <source point set> (+ options)         |\n"  );
  printf("  o------------------------------------------------------------------------o\n"  );
  printf("  ** Tab-separated files only. Extension of the input file MUST be '.txt'.  \n\n");
  printf(" OPTIONS:                                                                   \n"  );
  printf("  Parameters    -w <omega>, -l <lambda>, -g <gamma> -k <kappa>              \n"  );
  printf("  Acceleration  -J <rank:P>, -K <rank:G>, kdtree: -p, -r <seed>             \n"  );
  printf("  Downsampling  -D <type:x,y,b,X,Y,B, #points, radius>.                     \n"  );
  printf("  Kernel                                                                    \n"  );
  printf("    Standard:   -G <type:1,2,3> -b <beta>                                   \n"  );
  printf("    Geodesic:   -G'geo,<tau>,<triangles>'  | geo-kernel with a mesh         \n"  );
  printf("                -G'geo,<tau>,<nnk>,<nnr>'  | geo-kernel without a mesh      \n"  );
  printf("                -b <beta>, -z <epsilon>, -K <rank:G>                        \n"  );
  printf("  Convergence   -n <#loops:max>, -N <#loops:min>, -c <tolerance>            \n"  );
  printf("  Normalization -u <type;e,x,y,n>                                           \n"  );
  printf("  File Output   -o <prefix>, -s <variables:y,x,u,v,a,c,e,T,P,Y,A(=all)>     \n"  );
  printf("  Terminal I/O  quiet mode: -q, history mode: -h, warning-disabled mode: -W \n\n");
  printf("  *1) Parenthesis <...> specifies the argument of an option.                \n"  );
  printf("  *2) The downsampling option activates BCPD++/GBCPD++.                     \n\n");
  printf(" DEFAULT:                                                                   \n"  );
  printf("  -x X.txt, -y Y.txt, -w 0, -l 2, -b 2, -n 500, -o output_, -c 1e-4, -u e   \n"  );
  printf("  ** All accleration options are disabled unless specified explicitly.      \n\n");
  printf(" EXAMPLE:                                                                   \n"  );
  printf("  o------------------------------------------------------------------------o\n"  );
  printf("  | ./bcpd -x X.txt -y Y.txt -w0.1 -l2 -b2 -J300 -K80 -p -n90 -c1e-6 -svYP |\n"  );
  printf("  o------------------------------------------------------------------------o\n\n");
  printf(" REFERENCE:                                                                 \n"  );
  printf(" - Geodesic-Based Bayesian coherent point drift, IEEE TPAMI, 2022 (GBCPD).  \n"  );
  printf(" - A Bayesian formulation of coherent point drift, IEEE TPAMI, 2020 (BCPD). \n"  );
  printf(" - Accleelration of non-rigid point set registration with downsampling and  \n"  );
  printf("   Gaussian process regression, IEEE TPAMI, 2020 (BCDP++).                  \n\n");
  printf("\n");
}

static char warning[6][1024]= {
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  The numbers of points in both point sets are large. The execution |\n"
    "  |  will be considerably slow. Please use options, -J, -K, -p, and -D |\n"
    "  |  for accelerating the execution. For more detail, please read the  |\n"
    "  |  README.md file.                                                   |\n"
    "  o--------------------------------------------------------------------o\n\n"
    "  Do you want to continue the execution? [y/n] -- ",
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  The Nystrom option for computing P, i.e., '-J', is turned on, but |\n"
    "  |  '-p' is disabled. The computation will be quite inaccurate. The   |\n"
    "  |  use of the option '-p' in combination with '-J' is recommended    |\n"
    "  |  for more accurate computations.                                   |\n"
    "  o--------------------------------------------------------------------o\n",
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  The Nystrom option for computing P, i.e., '-J', is turned on, but |\n"
    "  |  the number of points J might be too small. The use of a larger J  |\n"
    "  |  is recommended for more accurate computations.                    |\n"
    "  o--------------------------------------------------------------------o\n",
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  The Nystrom option for computing G, i.e., '-K', is turned on, but |\n"
    "  |  the number of points K might be too small. The use of a larger K  |\n"
    "  |  is recommended for more accurate computations.                    |\n"
    "  o--------------------------------------------------------------------o\n",
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  You might want to solve regid registration for large point sets   |\n"
    "  |  because lambda is set to a very large value. If so, accelerate    |\n"
    "  |  carefully; use the following option:                              |\n"
    "  |                                                                    |\n"
    "  |    -J300 -K70 -p -d5 -e0.3 -f0.3 -g3 -D'B,2000,0.08' -w0.1,        |\n"
    "  |                                                                    |\n"
    "  |  for example. Otherwise, the computation will be unstable.         |\n"
    "  |  Do not output 'P', i.e., specify neither '-sP' nor '-sA' because  |\n"
    "  |  #nonzero elements in P will be enormous.                          |\n"
    "  o--------------------------------------------------------------------o\n",
    "  o--- !!! WARNING !!! ------------------------------------------------o\n"
    "  |  Excessively large beta might lead to an unstable computation due  |\n"
    "  |  to the rank deficiency of G. A smaller beta might be better.      |\n"
    "  o--------------------------------------------------------------------o\n",
};

static void yesno(void){
  char c; int n; n=scanf("%c",&c); if(n==0||(c!='y'&&c!='Y')){fprintf(stderr,"  Abort.\n\n"); exit(EXIT_SUCCESS);}
}

void printInfo(pwsz sz, pwpm pm){
  int nx,ny; double rx,ry; char stype[3][32]={"random sampling","inverse density","voxel grid"};
  nx=pm.dwn[TARGET]; rx=pm.dwr[TARGET];
  ny=pm.dwn[SOURCE]; ry=pm.dwr[SOURCE];

  fprintf(stderr,"\n  BCPD version %s (%s). OpenMP: %s.\n",_VERSION_,_DATE_,omp?"Turned on":"Not used");
  /* warning */
  if(pm.opt&PW_OPT_NWARN) goto skip;
  if(fmin(sz.M,ny)>=2000&&fmin(sz.N,nx)>=2000&&!sz.K) {fprintf(stderr,"\n%s",warning[0]);yesno();}
  if(pm.J&&!(pm.opt&PW_OPT_LOCAL))         {fprintf(stderr,"\n%s",warning[1]);}
  if(pm.J&&pm.J<150&&sz.M>=200&&sz.N>=200) {fprintf(stderr,"\n%s",warning[2]);}
  if(pm.K&&pm.K< 10&&sz.M>=100&&sz.N>=100) {fprintf(stderr,"\n%s",warning[3]);}
  if(pm.lmd>=1e6 &&sz.M>=2000&&sz.N>=2000) {fprintf(stderr,"\n%s",warning[4]);}
  if(pm.bet>=10&&pm.G<=3)                  {fprintf(stderr,"\n%s",warning[5]);}
  skip:

  fprintf(stderr,"\n");
  fprintf(stderr,"  Input Data:\n");
  fprintf(stderr,"    Point Set 1 (target): [%s]\n", pm.fn[TARGET]);
  fprintf(stderr,"    Point Set 2 (source): [%s]\n", pm.fn[SOURCE]);
  if(strlen(pm.fn[FACE_Y]))
  fprintf(stderr,"    Triangles   (source): [%s]\n", pm.fn[FACE_Y]);
  fprintf(stderr,"    Size of Point Set 1:  [%3d,%2d]\n", sz.N,sz.D);
  fprintf(stderr,"    Size of Point Set 2:  [%3d,%2d]\n", sz.M,sz.D);
  fprintf(stderr,"\n");

  fprintf(stderr,"  Parameters: \n");
  fprintf(stderr,"    omega   =  %.2lf\n", pm.omg);
  if(pm.lmd<1e3)   fprintf(stderr,"    lambda  =  %.2lf", pm.lmd);
  else             fprintf(stderr,"    lambda  =  %.2e",  pm.lmd);
  if(pm.G==0) fprintf(stderr,"  --> the expected drift length = %.3lf\n",sqrt(sz.D/pm.lmd));
  else        fprintf(stderr,"\n");
  fprintf(stderr,"    gamma   =  %.2lf\n", pm.gma);
  if(pm.kpa<=ZERO) fprintf(stderr,"    kappa   =  inf\n");
  else             fprintf(stderr,"    kappa   =  %lf\n", pm.kpa);

  fprintf(stderr,"\n");
  fprintf(stderr,"  Kernel:\n    ");
  switch(pm.G){
    case  0: if(pm.tau<=1e-5) fprintf(stderr,"[Gaussian]");
             else             fprintf(stderr,"[Geodesic]");
             break;
    case  1: fprintf(stderr,"[Inverse Multiquadric]"); break;
    case  2: fprintf(stderr,"[Rational Quadratic]");   break;
    case  3: fprintf(stderr,"[Laplacian]");            break;
    case  4: fprintf(stderr,"[Geodesic]");             break;
  }
  fprintf(stderr," with beta = %.2lf%s",pm.bet,pm.tau>1e-5?" and ":"\n");
  if(pm.tau>1e-5){ fprintf(stderr,"tau = %.2lf\n",pm.tau);
    if(pm.nnk>0) fprintf(stderr,"     *Surface Graph: %d-NNs, where max radius = %.3lf\n",pm.nnk,pm.nnr);
  }

  fprintf(stderr,"\n");
  fprintf(stderr,"  Acceleration: \n");
  if(nx||ny) fprintf(stderr,"    Downsampling:\n");
  if(nx){    fprintf(stderr,"      Target --> %5.d ",nx);
    if     (rx>0) fprintf(stderr,"[%s: r=%.3f]\n",stype[1], rx);
    else if(rx<0) fprintf(stderr,"[%s: r=%.3f]\n",stype[2],-rx);
    else          fprintf(stderr,"[%s]\n",stype[0]);
  }
  if(ny){    fprintf(stderr,"      Source --> %5.d ",ny);
    if     (ry>0) fprintf(stderr,"[%s: r=%.3f]\n",stype[1], ry);
    else if(ry<0) fprintf(stderr,"[%s: r=%.3f]\n",stype[2],-ry);
    else          fprintf(stderr,"[%s]\n",stype[0]);
  }
  if(sz.K)   fprintf(stderr,"    Fast G:  K = %d\n",sz.K);
  else       fprintf(stderr,"    Fast G:  OFF\n");
  if(sz.J||pm.opt&PW_OPT_LOCAL){
    fprintf(stderr,"    Fast P:  \n");
    if(sz.J) fprintf(stderr,"      Nystrom: J = %d\n", sz.J);
    else     fprintf(stderr,"      Nystrom: OFF\n");
    if(!(pm.opt&PW_OPT_LOCAL)) fprintf(stderr,"      KD tree: OFF\n");
    else fprintf(stderr,"      KD Tree: r = min(%.2lf,%.1lf*sigma) if sigma < %.2lf\n",pm.lim,pm.dlt,pm.btn);
  } else fprintf(stderr,"    Fast P:  OFF\n");
  if(pm.dwn[SOURCE]){
    if(pm.opt&PW_OPT_1NN) fprintf(stderr,"    Interpolation: 1NN\n");
    else {
      if(pm.K)  fprintf(stderr,"    Fast Interpolation:  L = %d\n",pm.K);
      else      fprintf(stderr,"    Fast Interpolation:  OFF\n");
    }
  }
  if(pm.rns&&(pm.J||pm.K)) fprintf(stderr,"    Rand Seed: %d\n",pm.rns);


  fprintf(stderr,"\n");
  fprintf(stderr,"  Convergence: \n");
  fprintf(stderr,"    tolerance   =  %.1e\n", pm.cnv);
  fprintf(stderr,"    max #loops  =  %d\n",   pm.nlp);
  fprintf(stderr,"    min #loops  =  %d\n",   pm.llp);

  fprintf(stderr,"\n");
  fprintf(stderr,"  Output: \n");
  fprintf(stderr,"    [%s*.txt]; * =",pm.fn[OUTPUT]);
  if(pm.opt&PW_OPT_SAVEY) fprintf(stderr," y");
  if(pm.opt&PW_OPT_SAVEX) fprintf(stderr,",x");
  if(pm.opt&PW_OPT_SAVEU) fprintf(stderr,",u");
  if(pm.opt&PW_OPT_SAVEV) fprintf(stderr,",v");
  if(pm.opt&PW_OPT_SAVEA) fprintf(stderr,",a");
  if(pm.opt&PW_OPT_SAVEC) fprintf(stderr,",c");
  if(pm.opt&PW_OPT_SAVEE) fprintf(stderr,",e");
  if(pm.opt&PW_OPT_SAVET) fprintf(stderr,",s,R,t");
  if(pm.opt&PW_OPT_SAVEP) fprintf(stderr,",P");
  if((pm.opt&PW_OPT_SAVES)&&(pm.opt&PW_OPT_DBIAS)) fprintf(stderr,",Sigma");
  fprintf(stderr,"\n\n");
}
