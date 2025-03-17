// Copyright (c) 2018-2024 Osamu Hirose
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

void gaussprod_l(
       double        *w,     /*  O  | M or N            | w=P1 or q                 */
       double        *U,     /*  O  | M x  D            | PX                        */
       double        *V,     /*  O  | M x  Df           | Pfx                       */
       double        *wd,    /*  W  |   *               | memory                    */
       int           *wi,    /*  W  |   *               | memory                    */
       const double  *Y,     /*  I  | D  x M            | input matrix Y            */
       const double  *X,     /*  I  | D  x N            | input matrix X            */
       const double  *fy,    /*  I  | Df x M            | function values fy        */
       const double  *fx,    /*  I  | Df x N            | function values fx        */
       const double  *q,     /*  I  | J = N or M        | weight vector q           */
       int           *T,     /* I/W | 3 x J +1          | kdtree                    */
       int            D,     /*  I  | const.            | dimension                 */
       int            Df,    /*  I  | const.            | dimension f               */
       int            M,     /*  I  | const.            | #points in Y              */
       int            N,     /*  I  | const.            | #points in X              */
       int            P,     /*  I  | const.            | #nystrom samples          */
       double         h,     /*  I  | const.            | gauss width for X, Y      */
       double        *hf,    /*  I  | const.            | gauss widths for fx, fy   */
       double         dlt,   /*  I  | const.            | neighbor width rate for h */
       double         lim,   /*  I  | const.            | maximum radius for kdtree */
       int            flg    /*  I  | const.            | flag:local+reuse+trans    */
  ){
  int d,p,i,j,k,I,J,L=M+N; const double *A,*B,*F,*G; int info; char uplo='U'; int *e=NULL;
  double *Z=NULL,*W=NULL,*K=NULL,*E=NULL; double reg=1e-8,rad=fmin(h*dlt,lim); int si=0,sd=0;
  int mtd=MAXTREEDEPTH; int *a,*u,*S,*bi; double *bd; double val;
  int tr=flg&GRAM_FLAG_TRANS; int nc=1+(tr?0:(D+Df));

  if(fx) assert(fy); else assert(fy==NULL&&Df==0);
  if(tr) assert(nc==1&&U==NULL&&V==NULL);

  /* memory: always use the same allocation if P>0 */
  if(P){e=wi+si;si+=L;Z=wd+sd;sd+=D*P;W=wd+sd;sd+=Df*P;K=wd+sd;sd+=P*P;E=wd+sd;sd+=P*nc;}
  /* init */
  if(tr) {I=N;J=M;A=X;B=Y;F=fx;G=fy;}
  else   {I=M;J=N;A=Y;B=X;F=fy;G=fx;}
  for(i=0;i<I;i++) w[i]=0;
  if(!tr){
    for(i=0;i<I;i++)for(d=0;d<D; d++) U[i+I*d]=0;
    for(i=0;i<I;i++)for(d=0;d<Df;d++) V[i+I*d]=0;
  }
  /* switch */
  if(flg&GRAM_FLAG_LOCAL) goto neighbor;
  if(P) goto nystrom; else goto direct;

  direct:
  #pragma omp parallel for private (j) private (d) private(val)
  for(i=0;i<I;i++)for(j=0;j<J;j++){
    val=q[j]*gauss(A+D*i,B+D*j,D,h)*(fx?anigauss(F+Df*i,G+Df*j,Df,hf):1);
    w[i]+=val; if(tr) continue;
    for(d=0;d<D; d++) U[i+I*d]+=val*B[d+D *j];
    for(d=0;d<Df;d++) V[i+I*d]+=val*G[d+Df*j];
  } return;

  neighbor: assert(T);
  bi=wi+si;wi+=6*J; a=wi+si;si+=I; S=wi+si;si+=mtd*I;
  bd=wd+sd;wd+=2*J; u=wi+si;si+=I;
  if(flg&GRAM_FLAG_BUILD) kdtree(T,bi,bd,B,D,J);
  #pragma omp parallel for private (j) private (d) private (val)
  for(i=0;i<I;i++){a[i]=u[i]=0;
    do{ eballsearch_next(a+i,S+mtd*i,u+i,A+D*i,rad,B,T,D,J); j=a[i];
      if(j>=0){ val=q[j]*gauss(A+D*i,B+D*j,D,h)*(fx?anigauss(F+Df*i,G+Df*j,Df,hf):1);
        w[i]+=val; if(tr) continue;
        for(d=0;d<D; d++) U[i+I*d]+=val*B[d+D *j];
        for(d=0;d<Df;d++) V[i+I*d]+=val*G[d+Df*j];
      }
    } while(u[i]);
  } return;

  nystrom: randperm(e,L);
  /* sampling */
  for(p=0;p<P;p++){k=e[p];
    if(k<M){/**/; for(d=0;d<D;d++){Z[d+D*p]=Y[d+D*k];} if(fy)for(d=0;d<Df;d++){W[d+Df*p]=fy[d+Df*k];}}
    else   {k-=M; for(d=0;d<D;d++){Z[d+D*p]=X[d+D*k];} if(fx)for(d=0;d<Df;d++){W[d+Df*p]=fx[d+Df*k];}}
  }
  /* E = KZB x q, qX, qfx */
  for(k=0;k<P*nc;k++) E[k]=0;
  #pragma omp parallel for private (j) private (d) private (val)
  for(p=0;p<P;p++)for(j=0;j<J;j++){
    val=q[j]*gauss(Z+D*p,B+D*j,D,h)*(fx?anigauss(W+Df*p,G+Df*j,Df,hf):1);
    E[p]+=val; if(tr) continue;
    for(d=0;d<D; d++) E[p+P*(d+1)  ]+=val*B[d+D *j];
    for(d=0;d<Df;d++) E[p+P*(d+1+D)]+=val*G[d+Df*j];
  }
  /* KZZ */
  for(i=0;i<P;i++)for(j=0;j<P;j++) K[i+P*j]=gauss(Z+D*i,Z+D*j,D,h)*(fx?anigauss(W+Df*i,W+Df*j,Df,hf):1)+(i==j?reg:0);
  /* inv(KZZ) x E */
  dpotrf_(&uplo,&P,K,&P,&info);          if(info!=0) goto err01;
  dpotrs_(&uplo,&P,&nc,K,&P,E,&P,&info); if(info!=0) goto err02;
  /* KAZ x E */
  #pragma omp parallel for private (p) private (d) private (val)
  for(i=0;i<I;i++)for(p=0;p<P;p++){
    val=gauss(A+D*i,Z+D*p,D,h)*(fy?anigauss(F+Df*i,W+Df*p,Df,hf):1);
    w[i]+=val*E[p]; if(tr) continue;
    for(d=0;d<D; d++) U[i+I*d]+=val*E[p+P*(d+1)  ];
    for(d=0;d<Df;d++) V[i+I*d]+=val*E[p+P*(d+1+D)];
  } 
  return;

  err01: printf("ERROR: The Cholesky factorization failed at gaussprod.\n"); exit(EXIT_FAILURE);
  err02: printf("ERROR: Solving linear equations failed at gaussprod.  \n"); exit(EXIT_FAILURE);
}

void gaussprod(
       double        *w,     /*  O  | M or N            | w=P1 or q                 */
       double        *U,     /*  O  | M x  D            | PX                        */
       double        *V,     /*  O  | M x  Df           | Pfx                       */
       double        *wd,    /*  W  |   *               | memory                    */
       int           *wi,    /*  W  |   *               | memory                    */
       const double  *Y,     /*  I  | D  x M            | input matrix Y            */
       const double  *X,     /*  I  | D  x N            | input matrix X            */
       const double  *fy,    /*  I  | Df x M            | function values fy        */
       const double  *fx,    /*  I  | Df x N            | function values fx        */
       const double  *q,     /*  I  | J = N or M        | weight vector q           */
       int           *T,     /* I/W | 3 x J +1          | kdtree                    */
       int            D,     /*  I  | const.            | dimension                 */
       int            Df,    /*  I  | const.            | dimension f               */
       int            M,     /*  I  | const.            | #points in Y              */
       int            N,     /*  I  | const.            | #points in X              */
       int            P,     /*  I  | const.            | #nystrom samples          */
       double         h,     /*  I  | const.            | gauss width for X, Y      */
       double        *hf,    /*  I  | const.            | gauss widths for fx, fy   */
       double         dlt,   /*  I  | const.            | neighbor width rate for h */
       double         lim,   /*  I  | const.            | maximum radius for kdtree */
       int            flg    /*  I  | const.            | flag:local+reuse+trans    */
  ){

  /*---------------+-----------------+--------------------+
  | storage:       |  int            |  double            |
  +----------------+-----------------+--------------------+
  | Nystrom        |  L              |  P^2+P(1+2Dc)      |
  | kd tree body   |  3L+1           |                    |
  | kd work build  |  6L             |  2L                |
  | kd work eball  |  (2+mtd)L       |                    |
  | whole space    |                 |   L x  3Dc         |
  +----------------+-----------------+-------------------*/

  int d,m,n,sd=0; int tr=flg&GRAM_FLAG_TRANS; double *lwd; double *x,*y,*Uc; int L=M+N,Dc=D+Df;

  // check arguments (for debugging)
  if(fy) assert(fx&&Df); else assert(fx==NULL&&Df==0);
  if(tr) assert(U==NULL&&V==NULL);

  /* case: truncated gauss (if local) */
  {gaussprod_l(w,U,V,wd,wi,Y,X,fy,fx,q,T,D,Df,M,N,P,h,hf,dlt,lim,flg); return;}

  /* case: barnes-hut kdtree (if local) */
  // allocation
  x=wd+sd;sd+=Dc*N; if(tr){Uc=NULL;} else {Uc=wd+sd;sd+=Dc*L;}
  y=wd+sd;sd+=Dc*M; lwd=wd+sd;
  // scaling (if Df>0 [function x whole] else [point set x whole])
  for(d=0;d<D; d++)for(n=0;n<N;n++) x[ d   +Dc*n]=X [d+D *n]/h;
  for(d=0;d<D; d++)for(m=0;m<M;m++) y[ d   +Dc*m]=Y [d+D *m]/h;
  for(d=0;d<Df;d++)for(n=0;n<N;n++) x[(d+D)+Dc*n]=fx[d+Df*n]/hf[d];
  for(d=0;d<Df;d++)for(m=0;m<M;m++) y[(d+D)+Dc*m]=fy[d+Df*m]/hf[d];

  // execution
  gaussprod_l(w,Uc,NULL,lwd,wi,y,x,NULL,NULL,q,T,Dc,0,M,N,P,1,NULL,dlt,lim,flg);
  // revert scaling
  if(tr) return;
  for(d=0;d<D; d++)for(m=0;m<M;m++) U[m+M*d]=Uc[m+M*d]*h;
  for(d=0;d<Df;d++)for(m=0;m<M;m++) V[m+M*d]=Uc[m+M*(d+D)]*hf[d];
}

