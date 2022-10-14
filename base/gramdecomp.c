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
#include<assert.h>
#include<math.h>
#include"misc.h"
#include"lapack.h"

void gramdecomp(
  double        * Q,                /* O |  M x K  | eigenvectors            */
  double        * L,                /* O |    K    | eigenvalues             */
  double        * wd,               /* W | K(K+11) | working memory (double) */
  int           * wi,               /* W |    M    | working memory (int)    */
  const double  * Y,                /* I |  D x M  | data matrix             */
  int             D,                /* I |    1    | dimension               */
  int             M,                /* I |    1    | #points                 */
  int             K,                /* I |    1    | #nystrom samples        */
  const double    bet,              /* I |    1    | kernel parameter        */
  double         (*kernel)(         /* I |   ptr   | kernel function         */
                   const double *,  /* I |    D    | point 1                 */
                   const double *,  /* I |    D    | point 2                 */
                   int,             /* I |    1    | D: dimension            */
                   double           /* I |    1    | kernel parameter        */
                 )
  ){

  int i,j,k,m; double c,val; double *A,*C,*Lr; int *U; int sd=0,si=0; char uplo='U',jobz='V'; int info,lwork; 
  assert(K>=1&&K<=M); 
  /* alias */
  lwork=K*10;
  A=wd+sd; sd+=K*K;   Lr=wd+sd; sd+=K;
  C=wd+sd; sd+=lwork; U =wi+si; si+=M;
  
  /* main computation */
  randperm(U,M); c=K/(double)M;
  for(i=0;i<K;i++)for(j=0;j<K;j++) A[i+K*j]=(*kernel)(Y+D*U[i],Y+D*U[j],D,bet);
  dsyev_(&jobz,&uplo,&K,A,&K,Lr,C,&lwork,&info); if(info!=0){goto err;}
  #pragma omp parallel for private (m) private (i) private(val)
  for(k=0;k<K;k++)for(m=0;m<M;m++)
    {val=0;for(i=0;i<K;i++){val+=A[i+K*(K-k-1)]*(*kernel)(Y+D*m,Y+D*U[i],D,bet);} Q[m+M*k]=sqrt(c)*val/Lr[K-k-1];}
  for(k=0;k<K;k++) L[k]=Lr[K-k-1]/c;
  for(k=0;k<K;k++) L[k]+=1e-9;

  return;
  err: printf("ERROR: The approximate eigendecomposition failed.\n"); exit(EXIT_FAILURE);

}
