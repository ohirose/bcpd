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

void gramdecomp(
  double        * Q,                /* O |  M x K  | eigenvectors            */
  double        * L,                /* O |    K    | eigenvalues             */
  double        * wd,               /* W | K(K+11) | working memory (double) */
  int           * wi,               /* W |    M    | working memory (int)    */
  const double  * Y,                /* I |  D x M  | data matrix             */
  int             D,                /* I |    1    | dimension               */
  int             M,                /* I |    1    | #points                 */
  int             K,                /* I |    1    | #nystrom samples        */
  const double    bet,              /* I |    1    | kernel parameter        */
  double         (*kernel)(         /* I |   ptr   | kernel function         */
                   const double *,  /* I |    D    | point 1                 */
                   const double *,  /* I |    D    | point 2                 */
                   int,             /* I |    1    | D: dimension            */
                   double           /* I |    1    | kernel parameter        */
                 )
);
