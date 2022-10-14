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
#include"heap.h"


#define QLOADED (*Q)

void dijkstra(
    double        * D,   /* O |    N    | distance          */
    int           * P,   /* O |    N    | previous node     */
    int           * wk,  /* W |  2N+1   | tag(N), heap(N+1) */
    const int    ** E,   /* I |  N x *  | edges             */
    const double ** W,   /* I |  N x *  | weights           */
    int             N,   /* I |  const. | #nodes            */
    int             s    /* I |  const. | start  node       */
  ){

  double d; int j,u,v; int *tag=wk,*Q=wk+N; assert(s<N);

  /* initialization */
  for(v=0;v<N;v++){D[v]=-1.0f;P[v]=-1;tag[v]=0;}
  heap_init(Q); 
  heap_insert(Q,D,s,0.0f);

  /* main computation */
  while(QLOADED){ heap_extract(&u,Q,D); tag[u]=1;
    for(j=1;j<=*(E[u]);j++){v=E[u][j]; if(tag[v]) continue; 
      d=D[u]+W[u][j]; /* temporal distance */
      if     (0>D[v]){P[v]=u;heap_insert (Q,D,v,d);} /* newly found v  */
      else if(d<D[v]){P[v]=u;heap_downkey(Q,D,v,d);} /* already found v */
    }
  }
}

static void swap   (int *a, int *b){ int tmp=*a; *a=*b; *b=tmp; }
static void reverse(int *a, int  N){ int n; for(n=0;n<N/2;n++) swap(a+n,a+N-1-n);}

void tracepath(int *path, int *len, const int *P, int N, int s, int g){
  int v=g,i=0;
  assert(g>=0&&g<N);
  assert(s>=0&&s<N);
  
  *len=0; for(*path=v;v!=-1;v=P[v]){path[++i]=P[v]; assert(v<N);}
  *len=i; assert(*len<=N);
  reverse(path,*len);
  assert(*path==s);
}

