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

void dijkstra(
    double        * D,   /* O |    N    | distance          */
    int           * P,   /* O |    N    | previous node     */
    int           * wk,  /* W |  2N+1   | tag(N), heap(N+1) */
    const int    ** E,   /* I |  N x *  | edges             */
    const double ** W,   /* I |  N x *  | weights           */
    int             N,   /* I |  const. | #nodes            */
    int             s    /* I |  const. | start  node       */
 );

void tracepath(int *path, int *len, const int *P, int N, int s, int g);
