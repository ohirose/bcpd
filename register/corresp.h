// Copyright (c) 2018-2026 Osamu Hirose
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

/* Include bcpd.h first. Returned arrays are freed by the caller. */
double *compute_hf(const double *e, pwsz sz, pwpm pm);
double *compute_lc(const double *X, const double *fx, const double *fx0, const double *e, double r, int N0, pwsz sz, pwpm pm);
double *extend_a(const double *Y0, const double *Y, const double *a, const int *U, int D, int M0, int M, int learned);
double *corresp_compute(
    const char   *prefix,
    const double *X,
    const double *y,
    const double *fx,
    const double *fy,
    const double *a,
    const double *lc,
    const double *hf,
    const double  r,
    pwsz          sz,
    pwpm          pm
);
