// Copyright (c) 2014--2019 Osamu Hirose
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

typedef unsigned char byte;

typedef struct sortbox   {double    val; int idx;} sortbox;

void prepare_sortbox      (sortbox *sb, const double * array, const int size);
int  cmp_sortbox_a        (const void *a, const void *b);
int  cmp_sortbox_d        (const void *a, const void *b);

char    ** calloc2c (const int uw, const int ulen);
double  ** calloc2d (const int M, const int N);
int     ** calloc2i (const int M, const int N);
double *** calloc3d (const int L, const int M, const int N);
short   ** calloc2s (const int M, const int N);
void       free2d   (double **a, int M);
void       free2i   (int    **a, int M);
void       fprint2d (const double **a, int M, int N);

double **  read2d   (int *nr, int *nc, char *mode, const char *file, const char *na);
void       write2d  (const char *file, const double **X, int nr, int nc, const char *fmt, const char *na);
void       conv2d   (const char *file, const char *fmt, const char *na);

char   ** read_strings (int *num, const char *file);
