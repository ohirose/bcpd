// Copyright (c) 2018-2019 Osamu Hirose
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

enum gram_flag{
  GRAM_FLAG_LOCAL = (1<<1),
  GRAM_FLAG_REUSE = (1<<2),
  GRAM_FLAG_TRANS = (1<<3),
  GRAM_FLAG_BUILD = (1<<4)
};

void gaussprod(
       double        *f,     /*  O  | I = M or N  | resulting values          */
       double        *wd,    /*  W  |   *         | memory                    */
       int           *wi,    /*  W  | M + N       | memory                    */
       const double  *Y,     /*  I  | D x M       | input matrix Y            */
       const double  *X,     /*  I  | D x N       | input matrix X            */
       const double  *q,     /*  I  | J = N or M  | weight vector q           */
       int           *T,     /* I/W | 3 x J +1    | kdtree                    */
       int            D,     /*  I  | const.      | dimension                 */
       int            M,     /*  I  | const.      | #points in Y              */
       int            N,     /*  I  | const.      | #points in X              */
       int            P,     /*  I  | const.      | #nystrom samples          */
       double         h,     /*  I  | const.      | band width                */
       double         dlt,   /*  I  | const.      | neighbor width rate for h */
       double         lim,   /*  I  | const.      | maximum radius for kdtree */
       int            flg    /*  I  | const.      | flag:local+reuse+trans    */
);

void gaussprod_batch(
       double        *w,     /*  O  | M or N            | required for w=P1         */
       double        *PX,    /*  O  | M x  D            | required for x=inv(w)     */
       double        *wd,    /*  W  |   *               | memory                    */
       int           *wi,    /*  W  |   *               | memory                    */
       const double  *Y,     /*  I  | D  x M            | input matrix Y            */
       const double  *X,     /*  I  | D  x N            | input matrix X            */
       const double  *q,     /*  I  | J = N or M        | weight vector q           */
       int           *T,     /* I/W | 3 x J +1          | kdtree                    */
       int            D,     /*  I  | const.            | dimension                 */
       int            M,     /*  I  | const.            | #points in Y              */
       int            N,     /*  I  | const.            | #points in X              */
       int            P,     /*  I  | const.            | #nystrom samples          */
       double         h,     /*  I  | const.            | gauss width for X, Y      */
       double         dlt,   /*  I  | const.            | neighbor width rate for h */
       double         lim,   /*  I  | const.            | maximum radius for kdtree */
       int            flg    /*  I  | const.            | flag:local+reuse+trans    */
);

void gaussprodf(
       double        *q,     /*  O  | M or N            | w=P1 or q                 */
       double        *wd,    /*  W  |   *               | memory                    */
       int           *wi,    /*  W  |   *               | memory                    */
       const double  *Y,     /*  I  | D  x M            | input matrix Y            */
       const double  *X,     /*  I  | D  x N            | input matrix X            */
       const double  *fy,    /*  I  | Df x M            | function values fy        */
       const double  *fx,    /*  I  | Df x N            | function values fx        */
       const double  *b,     /*  I  | J = N or M        | weight vector q           */
       int           *T,     /* I/W | 3 x J +1          | kdtree                    */
       int            D,     /*  I  | const.            | dimension                 */
       int            Df,    /*  I  | const.            | dimension f               */
       int            M,     /*  I  | const.            | #points in Y              */
       int            N,     /*  I  | const.            | #points in X              */
       int            P,     /*  I  | const.            | #nystrom samples          */
       double         h,     /*  I  | const.            | gauss width for X, Y      */
       double         hf,    /*  I  | const.            | gauss width for fx, fy    */
       double         dlt,   /*  I  | const.            | neighbor width rate for h */
       double         lim,   /*  I  | const.            | maximum radius for kdtree */
       int            flg    /*  I  | const.            | flag:local+reuse+trans    */
);

void gaussprodf_batch(
       double        *w,     /*  O  | M or N            | w                         */
       double        *PX,    /*  O  | M x  D            | PX                        */
       double        *Pfx,   /*  O  | M x  Df           | Pfx                       */
       double        *wd,    /*  W  |   *               | memory                    */
       int           *wi,    /*  W  |   *               | memory                    */
       const double  *Y,     /*  I  | D  x M            | input matrix Y            */
       const double  *X,     /*  I  | D  x N            | input matrix X            */
       const double  *fy,    /*  I  | Df x M            | function values fy        */
       const double  *fx,    /*  I  | Df x N            | function values fx        */
       const double  *q,     /*  I  | J = N or M        | weight vector q           */
       int           *T,     /* I/W | 3 x J +1          | kdtree                    */
       int            D,     /*  I  | const.            | dimension                 */
       int            Df,    /*  I  | const.            | dimension f               */
       int            M,     /*  I  | const.            | #points in Y              */
       int            N,     /*  I  | const.            | #points in X              */
       int            P,     /*  I  | const.            | #nystrom samples          */
       double         h,     /*  I  | const.            | gauss width for X, Y      */
       double         hf,    /*  I  | const.            | gauss width for fx, fy    */
       double         dlt,   /*  I  | const.            | neighbor width rate for h */
       double         lim,   /*  I  | const.            | maximum radius for kdtree */
       int            flg    /*  I  | const.            | flag:local+reuse+trans    */
);
