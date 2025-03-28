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
#include<stdio.h>
#define SQ(x) ((x)*(x))
double dist1    (const double *x, const double *y, int D){int d; double val=0;for(d=0;d<D;d++)val+=fabs(x[d]-y[d]); return val;}
double dist2    (const double *x, const double *y, int D){int d; double val=0;for(d=0;d<D;d++)val+=  SQ(x[d]-y[d]); return val;}
double sgauss   (const double *x, const double *y, int D){return exp(-dist2(x,y,D)/2.0);}
double gauss    (const double *x, const double *y, int D, double h){return exp(-dist2(x,y,D)/(2*SQ(h)));}
double laplace  (const double *x, const double *y, int D, double h){return exp(-dist1(x,y,D)/(h));}
double imquad   (const double *x, const double *y, int D, double h){return 1.0f/sqrt(SQ(h)+dist2(x,y,D));}
double rational (const double *x, const double *y, int D, double h){double val=dist2(x,y,D); return 1.0f-val/(val+SQ(h));}
double iwdist   (const double *x, const double *y, int D, double *w){int d; double val=0; for(d=0;d<D;d++){val+=SQ((x[d]-y[d])/w[d]);} return val;}
double anigauss (const double *x, const double *y, int D, double *w){return exp(-iwdist(x,y,D,w)/2);}
