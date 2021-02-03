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

#include<math.h>
#include<assert.h>

#define SQ(x) ((x)*(x))

double genrand64_real2 (void);

void   shuffle  (int *a, int n){int i,j,t; for(i=n-1;i;i--){j=(int)((i+1)*genrand64_real2()); t=a[i];a[i]=a[j];a[j]=t;}}
void   randperm (int *a, int n){int i; for(i=0;i<n;i++) a[i]=i; shuffle(a,n);}

double volume(const double *x, int D, int N){
  int d,n; double max,min,V=1.0;
  for(d=0;d<D;d++){
    max=x[d];for(n=1;n<N;n++) max=fmax(max,x[d+D*n]);
    min=x[d];for(n=1;n<N;n++) min=fmin(min,x[d+D*n]);
    V*=(max-min)*(N+1)/(double)N;
  }
  return V;
}

double det(const double *A, const int D){
  assert(D==2||D==3);
  return D==2? A[0]*A[3]-A[1]*A[2]:
               A[0+D*0]*(A[1+D*1]*A[2+D*2]-A[1+D*2]*A[2+D*1])
              -A[0+D*1]*(A[1+D*0]*A[2+D*2]-A[1+D*2]*A[2+D*0])
              +A[0+D*2]*(A[1+D*0]*A[2+D*1]-A[1+D*1]*A[2+D*0]);
}

