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
#include<assert.h>
#include<math.h>
#include"util.h"
#include"misc.h"
#include"kdtree.h"
#include"sgraph.h"
#include"lapack.h"


static sgraph* sgraph_new(int N){
  int n,ini=16;

  sgraph *sg=malloc(sizeof(sgraph));
  sg->N    = N;
  sg->E    = (int**)    calloc(N,sizeof(int*)); 
  sg->W    = (double**) calloc(N,sizeof(double*)); 
  sg->capa = (int*)     calloc(N,sizeof(int)); 

  for(n=0;n<N;n++){
    sg->E[n]=calloc(ini,sizeof(int));
    sg->W[n]=calloc(ini,sizeof(double));
    sg->capa[n]=ini;
  }

  return sg;
}

void sgraph_free(sgraph *sg){
  int n; int N=sg->N;  
  for(n=0;n<N;n++){free(sg->E[n]);free(sg->W[n]);}
  free(sg->capa); free(sg);
}


static int findedge(const int **E, int N, int from, int to){
  int j,nn=*(E[from]); 
  assert(from>=0&&from<N);
  assert(to  >=0&&to  <N);
  for(j=1;j<=nn&&E[from][j]!=to;j++){}
  return j>nn?0:j;
}

static int issymmetry(const int **E, const double **W, int N){
  int i,j,u,v; double w;
 
  for(u=0;u<N;u++)for(i=1;i<=*(E[u]);i++){ 
    v=E[u][i]; assert(v>=0&&v<N);
    w=W[u][i]; assert(w>=0);
    if(!(j=findedge(E,N,v,u))||!(w==W[v][j])) return 0; 
  }

  return 1;
}

/* add a directed edge; #nodes is assumed to be known */
static void add_edge(sgraph *sg, int from, int to, double w){ 
  int **E=sg->E,*capa=sg->capa,N=sg->N; double **W=sg->W; 
  int j,nn=*(E[from]),sz=capa[from]; int si=sizeof(int),sd=sizeof(double);

  assert(from>=0&&from<N);
  assert(to  >=0&&to  <N);
  assert(sz);

  /* add an edge only if it isn't included */
  j=findedge((const int**)E,N,from,to); if(j) return;

  /* reallocate an array if needed */
  if(1+nn==sz){sz*=2; capa[from]=sz;
    E[from]=realloc(E[from],sz*si);
    W[from]=realloc(W[from],sz*sd);
  } 

  nn++; 
  *(E[from])=nn; E[from][nn]=to;
  *(W[from])=nn; W[from][nn]=w;
}

/* add an undirected edge; #nodes is assumed to be known */
static void add_uedge(sgraph* sg, int u, int v, double w){ 
  add_edge(sg,u,v,w);
  add_edge(sg,v,u,w);
}

#define SQ(x) ((x)*(x))
static double dist(const double *x, const double *y, double D, int L){
  int d; double val=0;
  for(d=0;d<D;d++) val+=SQ(x[d]-y[d]);
  return L==1?sqrt(val):val;
}

sgraph* sgraph_from_points(const double *Y, int D, int M, int K, double emax){
  sgraph* sg; int j,u,v; int *T,**Q;

  Q =calloc2i(M,K+1);
  T =kdtree_build(Y,D,M);
  sg=sgraph_new(M);

  #pragma omp parallel for 
  for(v=0;v<M;v++) knnsearch(Q[v],K,emax,Y+D*v,v,Y,T,D,M);
  for(v=0;v<M;v++){K=*(Q[v]);
    for(j=1;j<=K;j++){u=Q[v][j];
      add_uedge(sg,v,u,dist(Y+D*v,Y+D*u,D,1));
    }
  } assert(issymmetry((const int**)sg->E,(const double**)sg->W,M));

  free(T);free2i(Q,M);
  return sg;
}

sgraph* sgraph_from_mesh(const double *Y, int D, int M, const char *file){
  sgraph *sg; int j,l,min; double **buff; int **line,nline,nc; char mode;

  /* read file */
  buff=read2d(&nline,&nc,&mode,file,"NA"); if(nc!=2&&nc!=3) goto err01;
  line=calloc2i(nline,nc); for(l=0;l<nline;l++)for(j=0;j<nc;j++) line[l][j]=(int)buff[l][j];
  free2d(buff,nline);

  /* check indices */
  min=line[0][0]; for(j=0;j<nc;j++)for(l=0;l<nline;l++){min=line[l][j]<min?line[l][j]:min;}
  if(min>1||min<0) goto err02;

  /* construct graph */
  sg=sgraph_new(M); sg->beg=min;
  for(l=0;l<nline;l++){ int v0,v1,v2;
    /* case: mesh & edge */
    v0=line[l][0]-(min>0?1:0); /* shift index */
    v1=line[l][1]-(min>0?1:0); /* shift index */
    add_uedge(sg,v0,v1,dist(Y+D*v0,Y+D*v1,D,1)); if(nc==2) continue;
    /* case: mesh */
    v2=line[l][2]-(min>0?1:0); /* shift index */
    add_uedge(sg,v1,v2,dist(Y+D*v1,Y+D*v2,D,1));
    add_uedge(sg,v2,v0,dist(Y+D*v2,Y+D*v0,D,1));
  } assert(issymmetry((const int**)sg->E,(const double**)sg->W,M));

  free2i(line,nline);
  return sg;

  err01: printf("ERROR: A edge/face definition file must contain 2 or 3 columns. Abort.\n"); exit(EXIT_FAILURE);
  err02: printf("ERROR: The minimum node index must be 0 or 1. Abort.\n");                   exit(EXIT_FAILURE);
}

