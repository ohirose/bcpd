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

#define ZERO 1e-8

typedef struct pairwise_size  {int D; int M; int N; int K; int J;} pwsz;
typedef struct pairwise_param {char fn[6][256]; char nrm; int opt;
  int nlp; int G; double dlt; double omg; double gma;
  int llp; int J; double lim; double lmd; double bet;
  int rns; int K; double btn; double kpa; double cnv;
  int dwn[2]; double dwr[2]; int nnk; double nnr; double tau; double eps;
} pwpm;

enum pairwise_option {
  PW_OPT_LOCAL = (1 << 0), PW_OPT_SAVEC = (1 << 9), PW_OPT_ACCEL = (1 <<18),
  PW_OPT_QUIET = (1 << 1), PW_OPT_SAVEE = (1 <<10), PW_OPT_PFLOG = (1 <<19),
  PW_OPT_DBIAS = (1 << 2), PW_OPT_SAVET = (1 <<11), PW_OPT_VTIME = (1 <<20),
  PW_OPT_HISTO = (1 << 3), PW_OPT_SAVEA = (1 <<12), PW_OPT_NWARN = (1 <<21),
  PW_OPT_SAVE  = (1 << 4), PW_OPT_SAVEP = (1 <<13), PW_OPT_1NN   = (1 <<22),
  PW_OPT_SAVEX = (1 << 5), PW_OPT_SAVES = (1 <<14), PW_OPT_INTPX = (1 <<23),
  PW_OPT_SAVEY = (1 << 6), PW_OPT_PATHX = (1 <<15), PW_OPT_NOSIM = (1 <<24),
  PW_OPT_SAVEU = (1 << 7), PW_OPT_PATHY = (1 <<16),
  PW_OPT_SAVEV = (1 << 8), PW_OPT_INFO  = (1 <<17),
};

enum {SOURCE=0, TARGET=1, OUTPUT=2, FACE_Y=3, FUNC_Y=4, FUNC_X=5, };

int bcpd(
  double        *   x,    /*  O  | DM x 1 (+nlp) | aligned target shape     */
  double        *   y,    /*  O  | DM x 1 (+nlp) | deformed source shape    */
  double        *   u,    /*  O  | DM x 1        | normalized def. shape    */
  double        *   v,    /*  O  | DM x 1        | displacement vectors     */
  double        *   w,    /*  O  |  M x 1        | #matches for each m      */
  double        *   a,    /*  O  |  M x 1        | mixing coefficients      */
  double        *   sgm,  /*  O  |  M x 1        | posterior covariance     */
  double        *   s,    /*  O  |    1          | scale factor             */
  double        *   R,    /*  O  |  D x D        | rotation matrix          */
  double        *   t,    /*  O  |    D          | translation vector       */
  double        *   r,    /*  O  |    1          | residual s.d.            */
  double        *   Np,   /*  O  |    1          | #matched points (est'd)  */
  double        *   pf,   /*  O  | nlp x 3       | comp. time (r/c) & sigma */
  double        *   wd,   /*  W  |    *          | working memory (double)  */
  int           *   wi,   /*  W  |    *          | working memory (int)     */
  const double  *   X,    /*  I  | DN x 1        | target point set         */
  const double  *   Y,    /*  I  | DM x 1        | source point set         */
  const double  *   LQ,   /*  I  | K + M x K     | only for geodesic kernel */
  const pwsz        sz,   /*  I  |               | D, M, N, K, J            */
  const pwpm        pm    /*  I  |               | tuning parameters        */
);

void interpolate(
  double         *T,
  const double   *Y,
  const int       N,
  const double   *x,
  const double   *y,
  const double   *w,
  const double   *s,
  const double   *R,
  const double   *t,
  const double   *r,
  const pwsz     sz,
  const pwpm     pm
);

void interpolate_1nn(
  double         *T,
  const double   *Y,
  const int       N,
  const double   *v,
  const double   *y,
  const double   *s,
  const double   *R,
  const double   *t,
  const pwsz     sz,
  const pwpm     pm
);

void interpolate_geok(
  double         *T,
  const double   *Y,
  const int       N,
  const double   *x,
  const double   *y,
  const double   *w,
  const double   *s,
  const double   *R,
  const double   *t,
  const double   *r,
  const double   *LQ,
  const int      *U,
  const pwsz     sz,
  const pwpm     pm
);

void interpolate_x(
    double       *x,
    const double *y,
    const double *X,
    int           D,
    int           M,
    int           N,
    const double  r,
    pwpm          pm
);
