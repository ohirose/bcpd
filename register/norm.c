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
#include<math.h>
#include<string.h>

#define SQ(x) ((x)*(x))

static void mean(double *mu, const double *X, int N, int D){
  int n,d; for(d=0;d<D;d++){mu[d]=0;for(n=0;n<N;n++){mu[d]+=X[d+D*n];} mu[d]/=N;}
}

static double scale(const double *X, const double *mu, int N, int D){
  int n,d; double val=0; for(d=0;d<D;d++)for(n=0;n<N;n++){val+=SQ(X[d+D*n]-mu[d]);}
  val/=N*D; return sqrt(val);
}

static void translate(double *X, const double *mu, int N, int D, int sign){
  int n,d; if(!mu) return;
  for(d=0;d<D;d++)for(n=0;n<N;n++) X[d+D*n]+=sign*mu[d];
}

static void resize(double *X, double sc, int N, int D, int inv){
  int n,d;
  if(inv) for(d=0;d<D;d++)for(n=0;n<N;n++){X[d+D*n]/=sc;}
  else    for(d=0;d<D;d++)for(n=0;n<N;n++){X[d+D*n]*=sc;}
}

static void norm_l(double *X, const double *mu, double sc, int N, int D, int revflag){
  int inv=1,plus=1,minus=-1;
  switch(revflag){
    case 0: translate(X,mu,N,D,minus); resize(X,sc,N,D,inv); break;
    case 1: resize(X,sc,N,D,0); translate(X,mu,N,D,plus);    break;
  }
}

/* alias */
void denormlize(double *X, const double *mu, double sc, int N, int D){norm_l(X,mu,sc,N,D,1);}

void normalize_batch(double *X, double *muX, double *scX, double *Y, double *muY, double *scY, int N, int M, int D, const char type){
  if(type=='n') return;
  mean(muX,X,N,D); *scX=scale(X,muX,N,D);
  mean(muY,Y,M,D); *scY=scale(Y,muY,M,D);

  if(type=='e'){norm_l(X,muX,*scX,N,D,0);norm_l(Y,muY,*scY,M,D,0);}
  if(type=='x'){norm_l(X,muX,*scX,N,D,0);norm_l(Y,muX,*scX,M,D,0);}
  if(type=='y'){norm_l(X,muY,*scY,N,D,0);norm_l(Y,muY,*scY,M,D,0);}
}

void denormalize_batch(double *X, const double *muX, double scX, double *Y, const double *muY, double scY, int N, int M, int D, const char type){
  int rev=1;
  if(type=='n') return;
  if(type=='e'){norm_l(X,muX,scX,N,D,rev);norm_l(Y,muX,scX,M,D,rev);}
  if(type=='x'){norm_l(X,muX,scX,N,D,rev);norm_l(Y,muX,scX,M,D,rev);}
  if(type=='y'){norm_l(X,muY,scY,N,D,rev);norm_l(Y,muY,scY,M,D,rev);}
}

void normalize(double *X, double *mu, double *sc, int N, int D){
  mean(mu,X,N,D);*sc=scale(X,mu,N,D);norm_l(X,mu,*sc,N,D,0);
}

