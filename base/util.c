// Copyright (c) 2014--2019 Osamu Hirose
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
#include<string.h>
#include<math.h>
#include<assert.h>
#include"util.h"

#define MAXC  100000
#define INIT  1024
#define MAXS  2048

double *** calloc3d (const int L, const int M, const int N){
  int l;
  double ***        a    = calloc   (L, sizeof(double**));
  for (l=0;l<L;l++) a[l] = calloc2d (M, N);
  return a;
}

void prepare_sortbox(sortbox *sb, const double * array, const int size){
  int i; for(i=0;i<size;i++){sb[i].val=array[i];sb[i].idx=i;}
  return;
}

int cmp_sortbox_a(const void *a, const void *b){
  const sortbox *sa = (const sortbox*) a;
  const sortbox *sb = (const sortbox*) b;
  return sa->val>sb->val ? 1 : sa->val<sb->val ? -1 : 0;
}

int cmp_sortbox_d(const void *a, const void *b){
  const sortbox *sa = (const sortbox*) a;
  const sortbox *sb = (const sortbox*) b;
  return sa->val<sb->val ? 1 : sa->val>sb->val ? -1 : 0;
}


double ** calloc2d (const int M, const int N){
  int m;
  double **         a    = calloc (M, sizeof(double*));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(double ));
  return a;
}

int    ** calloc2i (const int M, const int N){
  int m;
  int **            a    = calloc (M, sizeof(int *));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(int  ));
  return a;
}

short ** calloc2s (const int M, const int N){
  int m;
  short **          a    = calloc (M, sizeof(short*));
  for (m=0;m<M;m++) a[m] = calloc (N, sizeof(short ));
  return a;
}

char ** calloc2c (const int uw, const int ulen){
  int w;
  char **            a    = calloc (uw,  sizeof(char*));
  for (w=0;w<uw;w++) a[w] = calloc (ulen,sizeof(char ));
  return a;
}

void free2d(double **a, int M){int m;for(m=0;m<M;m++){free(a[m]);} free(a);}
void free2i(int    **a, int M){int m;for(m=0;m<M;m++){free(a[m]);} free(a);}

void write2d(const char *file, const double **X, int nr, int nc, const char* fmt, const char *na){
  FILE *fp; char *ext, mode='?'; double *buf; int m,n,M,N; double val; char str1[64],str2[64];
  size_t sz, si=sizeof(int),sd=sizeof(double); M=nr;N=nc;
  strcpy(str1,fmt); strcat(str1,"%c");
  strcpy(str2,na);  strcat(str2,"%c");

  if(!(ext=strrchr(file,'.'))){goto err00;} ext++;
  if(strcmp(ext,"bin")==0) mode='b';
  if(strcmp(ext,"txt")==0) mode='t';

  switch(mode){
    case 't': /* Tab-delimited file */
      fp=fopen(file,"w"); if(!fp) goto err01;
      for(m=0;m<M;m++)for(n=0;n<N;n++){val=X[m][n];
        if  (!isnan(val)) fprintf(fp,str1,X[m][n],n==N-1?'\n':'\t');
        else/*isnan(val)*/fprintf(fp,str2,        n==N-1?'\n':'\t');
      } fclose(fp); break;

    case 'b': /* Binary file */
      fp=fopen(file,"wb");if(!fp)  {goto err01;}
      if(1!=fwrite(&M,si,1,fp))    {goto err02;}
      if(1!=fwrite(&N,si,1,fp))    {goto err02;} sz=(size_t)N*M; buf=malloc(sd*sz);
      for(m=0;m<M;m++)for(n=0;n<N;n++) buf[n+N*m]=X[m][n];
      if(sz!=fwrite(buf,sd,sz,fp)) goto err02;
      free(buf); fclose(fp); break;

    case '?': goto err03;
  }

  return;

  err00: printf("ERROR: Missing extension: %s                  \n", file); exit(EXIT_FAILURE);
  err01: printf("ERROR: Failed to open: %s                     \n", file); exit(EXIT_FAILURE);
  err02: printf("ERROR: Failed to write: %s                    \n", file); exit(EXIT_FAILURE);
  err03: printf("ERROR: File '%s' may not be neither tab-",         file);
         printf("delimited file nor binary file. Abort.        \n");       exit(EXIT_FAILURE);
}

double ** read2d(int *nr, int *nc, char *mode, const char *file, const char *na){
  FILE *fp; char *s,*p,*ext,*dlm=",\t\n"; char *ptr;
  double *buf=NULL,**X; int m,n,M,N,lim=MAXC; size_t l,sz=INIT;
  size_t si=sizeof(int),sd=sizeof(double),sc=sizeof(char);

  *mode='?';
  if(!(ext=strrchr(file,'.'))){goto err00;} ext++;
  if(strcmp(ext,"bin")==0) *mode='b';
  if(strcmp(ext,"txt")==0) *mode='t';

  switch(*mode){
    case 't': /* Tab-delimited file */
      fp=fopen(file,"r"); if(!fp) goto err01;
      s  =malloc(sc*lim); /* storage line */
      buf=malloc(sd*sz);  /* storage data */
      l=0; /* initialize index */
      m=0; /* initialize #rows */
      /* for each line */
      while(fgets(s,lim,fp)){
        n=0; /* initialize #cols */
        m++; /* increment #rows  */
        if((1+strlen(s))==(size_t)lim) goto err02;
        /* for each token */
        for(p=strtok(s,dlm);p;p=strtok(NULL,dlm)){
          if   (0==strcmp(p,na))  {buf[l]=nan("NaN");}
          else/*0!=strcmp(p,na))*/{buf[l]=strtod(p,&ptr);if(p==ptr){goto err07;}}
          l++; /* increment index */
          n++; /* increment #cols */
          /* reallocation */
          if(l==sz){sz*=2;buf=realloc(buf,sd*sz);}
        } if(m==1) N=n; else if (N!=n) goto err03;
      } M=m; free(s); fclose(fp);break;

    case 'b': /* Binary file */
      fp=fopen(file,"rb");if(!fp) {goto err01;}
      if(1!=fread(&M,si,1,fp))    {goto err04;}
      if(1!=fread(&N,si,1,fp))    {goto err04;} sz=(size_t)N*M; buf=malloc(sd*sz);
      if(sz!=fread(buf,sd,sz,fp)) {goto err05;} fclose(fp);break;

    case '?': fp=fopen(file,"rb"); if(!fp) goto err01; else goto err06;
  } *nr=M;*nc=N; X=calloc2d(M,N);
  for(m=0;m<M;m++)for(n=0;n<N;n++) X[m][n]=buf[n+N*m];
  free(buf);

  return X;

  err00: printf("ERROR: Extension of the file '%s' is missing.  \n",    file); exit(EXIT_FAILURE);
  err01: printf("ERROR: File '%s' Not Found.                    \n",    file); exit(EXIT_FAILURE);
  err02: printf("ERROR: Number of characters in line %d of '%s' ",    m,file);
         printf("greater than the limit size. Abort.            \n");          exit(EXIT_FAILURE);
  err03: printf("ERROR: Number of columns are not identical in    ");
         printf("line %d of '%s'. Abort.                        \n",  m,file); exit(EXIT_FAILURE);
  err04: printf("ERROR: File '%s' may be empty. Abort.          \n",    file); exit(EXIT_FAILURE);
  err05: printf("ERROR: File '%s' may not be a matrix file.     \n",    file); exit(EXIT_FAILURE);
  err06: printf("ERROR: File '%s' may be neither a tab-delimited ",     file);
         printf("file nor a binary file. Abort.                 \n");          exit(EXIT_FAILURE);
  err07: printf("ERROR: Conversion from a string to a numeric "    );
         printf("value failed in line %d of '%s'. Abort.        \n",  m,file); exit(EXIT_FAILURE);
}

void conv2d(const char *file, const char *fmt, const char *na){
  int nr,nc; double **X; char mode,fout[256],*ext;

  X=read2d(&nr,&nc,&mode,file,na);

  strcpy(fout,file); ext=strrchr(fout,'.');
  if(mode=='t') strcpy(ext,".bin");
  if(mode=='b') strcpy(ext,".txt");

  write2d(fout,(const double**)X,nr,nc,fmt,na);
}

char** read_strings(int *num, const char *file){
  FILE *fp; char *s,**a; int l,sz,lim=MAXC,max=256;
  int sc=sizeof(char),ss=sizeof(char*);

  fp=fopen(file,"r"); if(!fp) goto err01;
  s =malloc(sc*lim);  /* storage: string       */
  a =malloc(ss*max);  /* storage: string array */
  l=0; /* initialize #lines */

  /* for each line */
  while(fgets(s,lim,fp)){
    sz=1+strlen(s); if(sz==lim) goto err02;
    /* skip if the first character is either # or % */
    if(strchr("#%\n",s[0])){continue;}
    a[l]=malloc(sc*sz); sscanf(s,"%s",a[l]);
    /* reallocation: array */
    if(l==max-1){max*=2;a=realloc(a,sc*max);}
    l++; /* increment #lines  */
  } *num=l;free(s);fclose(fp);

  return a;

  err01: printf("ERROR: Failed to open: %s\n",                        file);  exit(EXIT_FAILURE);
  err02: printf("ERROR: Number of characters in line %d of '%s' ",l+1,file);
         printf("greater than the limit size. Abort.\n"                   );  exit(EXIT_FAILURE);
}

