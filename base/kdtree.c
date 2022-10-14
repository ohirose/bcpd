// Copyright (c) 2017-2019 Osamu Hirose
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
#include"median.h"
#include"kdtree.h"

#define SQ(x) ((x)*(x))

static void swap(double *a, double *b){double tmp; tmp=*a; *a=*b; *b=tmp;}

static inline void divide(
   int          *  b,    /*  I/O | 5N | array to be divided         */
   double       *  v,    /*   W  | 2N | working memory              */
   int             K,    /*   I  |    | size of array to be divided */
   const double *  X,    /*   I  | DN | points                      */
   int             D,    /*   I  |    | dimension                   */
   int             N,    /*   I  |    | #points                     */
   int             p     /*   I  |    | current depth               */
  ){
  int k,i=0,j=0,e=0,c=K/2,d=p%D; double med,val; double *w=v+N;
  int *sz=b+N,*bl=b+2*N,*bc=b+3*N,*br=b+4*N;

  *sz=K; /* store original array size */
  for(k=0;k<K;k++){v[k]=X[d+D*b[k]];} swap(v,v+K/2);

  if (K==1) bl[0]=b[c];
  else { med=median(v,w,K);
    for(k=0;k<K;k++){ val=X[d+D*b[k]];
      if      (val<med) bl[i++]=b[k];
      else if (val>med) br[j++]=b[k];
      else              bc[e++]=b[k];
    }
    for(k=0;k<i;k++) b[k    ]=bl[k];
    for(k=0;k<e;k++) b[k+i  ]=bc[k];
    for(k=0;k<j;k++) b[k+i+e]=br[k];
  }
}

void kdtree(
       int           *T,  /* O | 3N+1 | depth(N),left(N),right(N),root(1)    */
       int           *a,  /* W |  6N  | index(N),size(N),buffer(3N),stack(N) */
       double        *v,  /* W |  2N  | buffer1(N), buffer2(N)               */
       const double  *X,  /* I |  DN  | points                               */
       int            D,  /* I |      | dimension                            */
       int            N   /* I |      | #points                              */
    ){
  int *bl=NULL,*br=NULL,*b=a,*S=a+5*N; int n,nl,nr,q=0,c,cl,cr,s,sl,sr; int p=0;

  /* init */
  for(n=0;n<3*N;n++) T[n]=-1;
  for(n=0;n<  N;n++) b[n]= n;
  /* basis */
  divide(b,v,N,X,D,N,p); S[q++]=0;c=N/2;T[b[c]]=p; T[3*N]=b[c]; /*root*/
  /* step */
  while(q){ b=a+S[--q];s=*(b+N); /*pop*/
    /* compute locations and sizes of divided arrays */
    c=s/2; n=b[c]; bl=b; br=bl+c+1; sl=c; sr=c-(s%2?0:1); p=T[n];
    /* division of coordinate space */
    if(sl){divide(bl,v,sl,X,D,N,p+1);S[q++]=bl-a;cl=sl/2;nl=bl[cl];T[nl]=p+1;T[n+N*1]=nl;}
    if(sr){divide(br,v,sr,X,D,N,p+1);S[q++]=br-a;cr=sr/2;nr=br[cr];T[nr]=p+1;T[n+N*2]=nr;}
  }
}

int *kdtree_build(const double *X, int D, int N){
  int *T,*wi; double *wd;

  T =calloc(3*N+1,sizeof(int));
  wi=calloc(6*N,  sizeof(int));
  wd=calloc(2*N,  sizeof(double));
  kdtree(T,wi,wd,X,D,N);
  free(wi); free(wd);

  return T;
}

static double dist(const double *x, const double *y, double D){
  int d; double val=0;
  for(d=0;d<D;d++) val+=SQ(x[d]-y[d]);
  return sqrt(val);
}

void eballsearch_next(
       int           *m,  /*  O  |   1   | next index within radius e */
       int           *S,  /*  W  | 2logN | stack                      */
       int           *q,  /*  W  |   1   | top index of stack 'S'     */
       const double  *y,  /*  I  |   D   | the point of interest      */
       double         e,  /*  I  | const.| ball radius                */
       const double  *X,  /*  I  |   DN  | points                     */
       const int     *T,  /*  I  |  3N+1 | kdtree                     */
       int            D,  /*  I  | const.| dimension                  */
       int            N   /*  I  | const.| #points                    */
     ){

  int d,p,nl,nr,n=T[3*N],state=1; double u,v;

  if(*q==0) S[(*q)++]=n;
  while((*q)&&state){ n=S[--(*q)];
     nl=T[n+N*1]; p=T[n];
     nr=T[n+N*2]; d=p%D;

     if(dist(y,X+D*n,D)<=e){*m=n;state=0;}

     v=y[d]-X[d+D*n]; u=fabs(v);
     if   (v>0){if(nr>=0)S[(*q)++]=nr; if(nl>=0&&u<=e)S[(*q)++]=nl;}
     else      {if(nl>=0)S[(*q)++]=nl; if(nr>=0&&u<=e)S[(*q)++]=nr;}

  }  if( state) *m=-1;
}

static void nn_approx(
       int           *i,  /*  O  | const.| subspace id               */
       double        *e,  /*  O  | const.| distance                  */
       const double  *y,  /*  I  | const.| the point of interest     */
       const double  *X,  /*  I  |  DxN  | points                    */
       const int     *T,  /*  I  | 3xN+1 | kdtree                    */
       int            D,  /*  I  | const.| dimension                 */
       int            N   /*  I  | const.| #points                   */
     ){

  int d,n=T[3*N],nl,nr; double v,l; *e=1e250;
  while(1){
     nl=T[n+N*1]; d=T[n]%D;
     nr=T[n+N*2]; v=y[d]-X[d+D*n];
     l=dist(y,X+D*n,D); if(l<*e){*e=l;*i=n;}
     if      (v> 0&&nr>=0) n=nr;
     else if (v<=0&&nl>=0) n=nl;
     else    break;
  }
}

void nnsearch(
       int           *i,  /*  O  | const.| nearest neighbor          */
       double        *e,  /*  O  | const.| distance                  */
       const double  *y,  /*  I  | const.| the point of interest     */
       const double  *X,  /*  I  |  DxN  | points                    */
       const int     *T,  /*  I  | 3xN+1 | kdtree                    */
       int            D,  /*  I  | const.| dimension                 */
       int            N   /*  I  | const.| #points                   */
     ){

  int d,p,q=0,n=T[3*N],nl,nr; int S[MAXTREEDEPTH]; double u,v,dst;

  nn_approx(i,e,y,X,T,D,N); S[q++]=n;
  while(q){assert(q>=0||(q&&S[q-1]>=0)); n=S[--q];
     nl=T[n+N*1]; p=T[n]; assert(p<MAXTREEDEPTH);
     nr=T[n+N*2]; d=p%D;

     dst=dist(y,X+D*n,D); if(dst<=*e){*e=dst;*i=n;}

     v=y[d]-X[d+D*n]; u=fabs(v);
     if   (v>0){if(nr>=0)S[q++]=nr; if(nl>=0&&u<=*e)S[q++]=nl;}
     else      {if(nl>=0)S[q++]=nl; if(nr>=0&&u<=*e)S[q++]=nr;}
  }
}

void knnsearch(
       int           *Q,  /*  O  |  1+K  | k nearest neighbors (+size) */
       int            K,  /* I/O | const.| #neighbors                  */
       double         e,  /* I/O | const.| maximum distance            */
       const double  *y,  /*  I  | const.| the point of interest       */
       int            me, /*  I  | const.| y's index if y in X else <0 */
       const double  *X,  /*  I  |  DxN  | points                      */
       const int     *T,  /*  I  | 3xN+1 | kdtree                      */
       int            D,  /*  I  | const.| dimension                   */
       int            N   /*  I  | const.| #points                     */
     ){

  int d,p,q=0,k=0,n=T[3*N],nl,nr,kmax=0,num; int S[MAXTREEDEPTH];
  double u,v,dst,val;

  num=0; for(k=1;k<K;k++) Q[k]=-1;

  S[q++]=n;
  while(q){assert(q>=0||(q&&S[q-1]>=0)); n=S[--q];
    nl=T[n+N*1]; p=T[n]; assert(p<MAXTREEDEPTH);
    nr=T[n+N*2]; d=p%D;

    dst=dist(y,X+D*n,D);
    if(dst<e&&n!=me){ /* exclude y if 'me' is non-negative */
      if(num< K){Q[++num]=n;}
      else      {Q[kmax ]=n;assert(kmax>0);} /* replace with current max */
      if(num==K) for(k=1;k<=K;k++){val=dist(y,X+D*Q[k],D);if(k==1||val>e){e=val;kmax=k;}} /* update max */
    }

    v=y[d]-X[d+D*n]; u=fabs(v);
    if   (v>0){if(nl>=0&&u<=e)S[q++]=nl; if(nr>=0)S[q++]=nr;}
    else      {if(nr>=0&&u<=e)S[q++]=nr; if(nl>=0)S[q++]=nl;}
  } *Q=num;
}

