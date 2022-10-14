// Copyright (c) 2021-2022 Osamu Hirose
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
#include"lapack.h"
#include"misc.h"
#include"sampling.h"
#include"dijkstra.h"

static double gaussd(double d, double h){return exp(-d*d/(2*h*h));}
static double dist (const double *x, const double *y, int D){int d; double val=0;for(d=0;d<D;d++)val+=(x[d]-y[d])*(x[d]-y[d]); return sqrt(val);}

double *geokdecomp(int *Knew, const double *Y, int D, int M, const int **E, const double **W, int K, double bet, double tau, double eps){
  int i,j,k,m; double *kerns; int *prevs,*works,*U; char uplo='U',jobz='V'; int info,lwork; double *dummy;
  double *A,*L,*Lr,*Q,*C,*LQ; double c,val; int si=sizeof(int),sd=sizeof(double); double a=tau;

  assert(K);

  lwork=K*10;
  kerns=calloc((  M  )*K,sd); L =calloc(M,sd); Q=calloc(M*K,  sd);
  prevs=calloc((  M  )*K,si); Lr=calloc(K,sd); A=calloc(K*K,  sd);
  works=calloc((2*M+1)*K,si); U =calloc(M,si); C=calloc(lwork,sd);

  /* geodesic kernel (partial) computation */
  dummy=calloc(D*K,sd); assert(dummy); downsample(dummy,U,K,(double*)Y,D,M,-0.05); free(dummy);
  #pragma omp parallel for
  for(j=0;j<K;j++) dijkstra(kerns+M*j,prevs+M*j,works+(2*M+1)*j,(const int**)E,(const double**)W,M,U[j]);
  for(j=0;j<K;j++)for(m=0;m<M;m++) kerns[m+M*j]=a*((kerns[m+M*j]<0)?0:gaussd(kerns[m+M*j],bet))+(1-a)*gaussd(dist(Y+D*m,Y+D*U[j],D),bet);

  /* approximate eigendecomposition */
  c=K/(double)M;
  for(i=0;i<K;i++)for(j=0;j<K;j++) A[i+K*j]=kerns[U[i]+M*j];
  dsyev_(&jobz,&uplo,&K,A,&K,Lr,C,&lwork,&info); if(info!=0){goto err;}

  #pragma omp parallel for private (m) private (i) private(val)
  for(k=0;k<K;k++)for(m=0;m<M;m++)
    {val=0;for(i=0;i<K;i++){val+=A[i+K*(K-k-1)]*kerns[m+M*i];} Q[m+M*k]=sqrt(c)*val/Lr[K-k-1];}
  for(k=0;k<K;k++) L[k]=Lr[K-k-1]/c;
  for(k=1;k<K&&L[k]/L[0]>eps;k++){} *Knew=k;

  /* copy to LQ */
  K=*Knew; LQ=calloc(K+M*K,sd);
  for(k=0;k<K;k++) LQ[k]=L[k];
  for(m=0;m<M;m++)for(k=0;k<K;k++) LQ[m+M*k+K]=Q[m+M*k];

  free(kerns); free(L);  free(Q);
  free(works); free(Lr); free(C);
  free(prevs); free(U);  free(A);

  return LQ;

  err: fprintf(stderr,"ERROR: dsyev in geokdecomp.c.\n"); exit(EXIT_FAILURE);
}

