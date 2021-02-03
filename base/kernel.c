// Copyright (c) 2019 Osamu Hirose
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
#define SQ(x) ((x)*(x))
double dist1    (const double *x, const double *y, int D){int d; double val=0;for(d=0;d<D;d++)val+=fabs(x[d]-y[d]); return val;}
double dist2    (const double *x, const double *y, int D){int d; double val=0;for(d=0;d<D;d++)val+=  SQ(x[d]-y[d]); return val;}
double gauss    (const double *x, const double *y, int D, double h){return exp(-dist2(x,y,D)/(2*h*h));}
double gaussk   (const double *x, const double *y, int D, const double *h){return exp(-dist2(x,y,D)/(2*SQ(*h)));}
double laplace  (const double *x, const double *y, int D, const double *h){return exp(-dist1(x,y,D)/(*h));}
double imquad   (const double *x, const double *y, int D, const double *h){return 1.0f/sqrt(SQ(*h)+dist2(x,y,D));}
double rational (const double *x, const double *y, int D, const double *h){double val=dist2(x,y,D); return 1.0f-val/(val+SQ(*h));}
double neural   (const double *x, const double *y, int D, const double *h){
  int d; double xy,xx,yy,val,c0=SQ(h[0]),c1=SQ(h[1]);
  val=c0;for(d=0;d<D;d++){val+=c1*x[d]*x[d];} xx=val;
  val=c0;for(d=0;d<D;d++){val+=c1*y[d]*y[d];} yy=val;
  val=c0;for(d=0;d<D;d++){val+=c1*x[d]*y[d];} xy=val;
  return (2/M_PI)*asin(2*xy/sqrt((1+2*xx)*(1+2*yy)));
}

/* You can add your own kernel by inserting its definition into the following code. */
/* At most two parameters of a kernel can be stored in *h.                          */
double mykernel  (const double *x, const double *y, int D, const double *h){
  /* BEGIN: This is a dummy code. Modify it as you like. */
  return exp(-dist2(x,y,D)/(2*SQ(*h)));
  /* END:   This is a dummy code. Modify it as you like. */
}
