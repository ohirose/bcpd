// Copyright (c) 2014 Osamu Hirose
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

double median(double *a, double *w, const int N){
  const int c = N/2;
  int     i,j,k,l,u,e,ofs=0,size=N;
  double  *tmp, p;/*pivot*/

  while(1){i=j=k=0;p=a[0];e=1;
    for(i=1;i<size;i++){
      if      (a[i]< p)  a[j++]=a[i];
      else if (a[i]> p)  w[k++]=a[i];
      else   /*a[i]==p*/ e++;
    } l=ofs+j;u=l+e-1;

    if      (c<l) {size=j;}
    else if (c>u) {tmp=a;a=w;w=tmp;ofs=u+1;size=k;}
    else break;
  }

  return p;
}
