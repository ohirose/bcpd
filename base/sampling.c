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
#include"kdtree.h"
#include"misc.h"
#define SQ(x) ((x)*(x))

double genrand64_real1(void);

static void resampling(int *nums, const double *probs, int N, int M){
  int i; double u =genrand64_real1()/(double)M;
  for(i=0;i<N;i++){
    nums[i]=floor((probs[i]-u)*M)+1; 
    u+=(nums[i]/(double)M)-probs[i];
  }
}

static void downsample_a(double *x, int *U, int L, double *X, int D, int N){
  int d,l; if(L>N){goto err01;}

  randperm(U,N); for(d=0;d<D;d++)for(l=0;l<L;l++){x[d+D*l]=X[d+D*U[l]];}
  return;

  err01:
  printf("\n\n  ERROR: L<=N must be satisfied in the function 'downsample_a'. Abort.\n\n");
  exit(EXIT_FAILURE);
}

static void downsample_b(double *x, int *U, int L, double *X, int D, int N, double e){
  int *S,*T,*a,*c; double *w,*v; int sd=sizeof(double),si=sizeof(int);
  int d,j,n,q,l=0,mtd=MAXTREEDEPTH; double val=0;
  /* allocation */
  T=calloc(3*N+1,si); S=calloc(N*mtd,si);
  a=calloc(6*N,  si); w=calloc(N,sd);
  v=calloc(2*N,  sd); c=calloc(N,si);
  /* build kdtree */ 
  kdtree(T,a,v,X,D,N);
  /* count #neighbors */
  #pragma omp parallel for private (j) private (q)
  for(n=0;n<N;n++){j=q=0;w[n]=0.0f;
    do{eballsearch_next(&j,S+mtd*n,&q,X+D*n,e,X,T,D,N);if(j>=0){w[n]+=1.0f;}} while(q);
    assert(w[n]>=1.0f);
  }
  /* sampling probabilities */
  for(n=0;n<N;n++) val+=1.0f/(w[n]);
  for(n=0;n<N;n++) w[n]=1.0f/(w[n]*val);
  /* resampling */
  resampling(c,w,N,L);
  /* output */
  for(n=0;n<N;n++)for(j=0;j<c[n];j++){U[l]=n;for(d=0;d<D;d++){x[d+D*l]=X[d+D*n];} l++;}

  free(T);free(a);free(v);
  free(S);free(w);free(c);
}

/* voxel grid filter */
static void downsample_c(double *x, int *U, int L, double *X, int D, int N, double e){
  int d,j,l=0,n,num; size_t K; int *v,*c,*np,*cum,*div; double *w,*max,*min; int sd=sizeof(double),si=sizeof(int);
  double val=0;
  /* allocation */
  v=calloc(N,si); max=calloc(D,sd); div=calloc(D,si); w=calloc(N,sd);
  c=calloc(N,si); min=calloc(D,sd); cum=calloc(D,si);
  /* bounding box */
  for(d=0;d<D;d++){min[d]=X[d];for(n=0;n<N;n++){min[d]=X[d+D*n]<min[d]?X[d+D*n]:min[d];}}
  for(d=0;d<D;d++){max[d]=X[d];for(n=0;n<N;n++){max[d]=X[d+D*n]>max[d]?X[d+D*n]:max[d];}}
  /* divide in grid & count points in a voxel */
  for(d=0;d<D;d++) div[d]=ceil((max[d]-min[d])/e);
  cum[0]=1; cum[1]=div[0]; for(d=2;d<D;d++) cum[d]=cum[d-1]*div[d-1];
  K=cum[D-1]*div[D-1]; if(K>=1e8){printf("  ERROR: Voxel grid width is too small. Abort.\n\n"); exit(EXIT_FAILURE);}
  np=calloc(K,si);
  for(n=0;n<N;n++){v[n]=0;for(d=0;d<D;d++){j=floor((X[d+D*n]-min[d])/e);j-=(j==div[d])?1:0;v[n]+=cum[d]*j;}}
  for(n=0;n<N;n++) np[v[n]]++;
  /* sampling probabilities */
  for(n=0;n<N;n++){num=np[v[n]];assert(num>0);w[n]=1.0f/num;val+=w[n];}
  for(n=0;n<N;n++) w[n]/=val;
  /* resampling */
  resampling(c,w,N,L);
  /* output */
  for(n=0;n<N;n++)for(j=0;j<c[n];j++){U[l]=n;for(d=0;d<D;d++){x[d+D*l]=X[d+D*n];} l++;}

  free(v);free(max);free(div);free(w);
  free(c);free(min);free(cum);free(np);
}

void downsample(double *x, int *U, int L, double *X, int D, int N, double e){
  if     (e<0) downsample_c(x,U,L,X,D,N,-e);
  else if(e>0) downsample_b(x,U,L,X,D,N, e);
  else         downsample_a(x,U,L,X,D,N);
}

