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

#define MAXTREEDEPTH 32

void kdtree(
       int           *T,  /* O | 3N+1 | depth(N),left(N),right(N),root(1)    */
       int           *a,  /* W |  6N  | index(N),size(N),stack(N),buffer(3N) */
       double        *v,  /* W |  2N  | buffer1(N), buffer2(N)               */
       const double  *X,  /* I |  DN  | points                               */
       int            D,  /* I |      | dimension                            */
       int            N   /* I |      | #points                              */
);

int *kdtree_build (const double *X, int D, int N);

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
);

void nnsearch(
       int           *i,  /*  O  | const.| nearest neighbor          */
       double        *e,  /*  O  | const.| distance                  */
       const double  *y,  /*  I  | const.| the point of interest     */
       const double  *X,  /*  I  |  DxN  | points                    */
       const int     *T,  /*  I  | 3xN+1 | kdtree                    */
       int            D,  /*  I  | const.| dimension                 */
       int            N   /*  I  | const.| #points                   */
);

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
);
