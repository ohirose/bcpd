// Copyright (c) 2018-2019 Osamu Hirose
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
#include<math.h>
#include<assert.h>
#include<time.h>
#include"misc.h"
#include"kdtree.h"
#include"lapack.h"
#include"kernel.h"
#include"gaussprod.h"

void gaussprod(
       double        *f,     /*  O  | I = M or N  | resulting values          */
       double        *wd,    /*  W  |   *         | memory                    */
       int           *wi,    /*  W  |   *         | memory                    */
       const double  *Y,     /*  I  | D x M       | input matrix Y            */
       const double  *X,     /*  I  | D x N       | input matrix X            */
       const double  *q,     /*  I  | J = N or M  | weight vector q           */
       int           *T,     /* I/W | 3 x J +1    | kdtree                    */
       int            D,     /*  I  | const.      | dimension                 */
       int            M,     /*  I  | const.      | #points in Y              */
       int            N,     /*  I  | const.      | #points in X              */
       int            P,     /*  I  | const.      | #nystrom samples          */
       double         h,     /*  I  | const.      | band width                */
       double         dlt,   /*  I  | const.      | neighbor width rate for h */
       double         lim,   /*  I  | const.      | maximum radius for kdtree */
       int            flg    /*  I  | const.      | flag:local+reuse+trans    */
  ){
  int d,m,n,p,i,j,I,J,L=M+N; const double *A,*B; int info,one=1; char uplo='U'; int *U=NULL;
  double *Z=NULL,*K=NULL,*v=NULL; double reg=1e-8,rad=fmin(h*dlt,lim); int si=0,sd=0;
  int mtd=MAXTREEDEPTH; int *a,*u,*S,*bi; double *bd;

  /* memory: always use the same allocation if P>0 */
  if(P){U=wi+si;si+=L; Z=wd+sd;sd+=D*L; v=wd+sd;sd+=P; K=wd+sd;sd+=P*P;}
  /* switch */
  if(flg&GRAM_FLAG_TRANS) {I=N;J=M;A=X;B=Y;} else   {I=M;J=N;A=Y;B=X;}
  if(flg&GRAM_FLAG_LOCAL) {goto neighbor;  } if(!P) {goto direct;    }
  if(flg&GRAM_FLAG_REUSE) {goto nystrom_r; } else   {goto nystrom;   }

  direct:
  #pragma omp parallel for private (j)
  for(i=0;i<I;i++){f[i]=0;for(j=0;j<J;j++)f[i]+=q[j]*gauss(A+D*i,B+D*j,D,h);} return;

  neighbor: assert(T);
  bi=wi+si;wi+=6*J; a=wi+si;si+=I; S=wi+si;si+=mtd*I;
  bd=wd+sd;wd+=2*J; u=wi+si;si+=I;
  if(flg&GRAM_FLAG_BUILD) kdtree(T,bi,bd,B,D,J);
  #pragma omp parallel for private (j)
  for(i=0;i<I;i++){a[i]=u[i]=f[i]=0;
    do{
      eballsearch_next(a+i,S+mtd*i,u+i,A+D*i,rad,B,T,D,J);
      j=a[i]; if(j>=0) f[i]+=q[j]*gauss(A+D*i,B+D*j,D,h);
    } while(u[i]);
  }
  return;

  nystrom: randperm(U,L);
  for(d=0;d<D;d++)for(m=0;m<M;m++) Z[d+D*( m )]=Y[d+D*m];
  for(d=0;d<D;d++)for(n=0;n<N;n++) Z[d+D*(M+n)]=X[d+D*n];
  for(i=0;i<P;i++)for(j=0;j<P;j++) K[i+P*j]=gauss(Z+D*U[i],Z+D*U[j],D,h)+(i==j?reg:0);
  dpotrf_(&uplo,&P,K,&P,&info); if(info!=0){goto err01;}

  nystrom_r:
  #pragma omp parallel for private (j)
  for(p=0;p<P;p++){v[p]=0;for(j=0;j<J;j++)v[p]+=gauss(Z+D*U[p],B+D*j,D,h)*q[j];}
  dpotrs_(&uplo,&P,&one,K,&P,v,&P,&info); if(info!=0){goto err02;}
  #pragma omp parallel for private (p)
  for(i=0;i<I;i++){f[i]=0;for(p=0;p<P;p++)f[i]+=gauss(Z+D*U[p],A+D*i,D,h)*v[p];}
  return;

  err01: printf("ERROR: The Cholesky factorization failed at gaussprod.\n"); exit(EXIT_FAILURE);
  err02: printf("ERROR: Solving linear equations failed at gaussprod.  \n"); exit(EXIT_FAILURE);
}

/* faster implementation with kdtree, specialized for w and PX */
void gaussprod_kdbatch(
       double        *w,     /*  O  |   M      | required for w=P1         */
       double        *PX,    /*  O  |   M      | required for x=inv(w)*PX  */
       int           *wi,    /*  W  |   *      | memory                    */
       const double  *Y,     /*  I  | D x M    | input matrix Y            */
       const double  *X,     /*  I  | D x N    | input matrix X            */
       const double  *q,     /*  I  |   N      | weight vector q           */
       int           *T,     /* I/W | 3 x N +1 | kdtree                    */
       int            D,     /*  I  | const.   | dimension                 */
       int            M,     /*  I  | const.   | #points in Y              */
       int            N,     /*  I  | const.   | #points in X              */
       double         h,     /*  I  | const.   | band width                */
       double         dlt,   /*  I  | const.   | neighbor width rate for h */
       double         lim    /*  I  | const.   | maximum radius for kdtree */
  ){
  int d,m,n; double rad=fmin(h*dlt,lim); int si=0; int mtd=MAXTREEDEPTH;
  double val; int *a,*u,*S;

  assert(T); a=wi+si;si+=M; u=wi+si;si+=M; S=wi+si;si+=mtd*M;
  #pragma omp parallel for private (d) private (n) private(val)
  for(m=0;m<M;m++){a[m]=u[m]=w[m]=0;for(d=0;d<D;d++)PX[m+M*d]=0;
    do{
      eballsearch_next(a+m,S+mtd*m,u+m,Y+D*m,rad,X,T,D,N); n=a[m];
      if(n>=0){val=q[n]*gauss(Y+D*m,X+D*n,D,h);w[m]+=val;for(d=0;d<D;d++)PX[m+M*d]+=X[d+D*n]*val;}
    } while(u[m]);
  }
  return;

}
