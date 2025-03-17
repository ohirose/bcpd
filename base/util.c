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

int charcount(FILE *fp){
  int c; int ct=0;
  fseek(fp,0,SEEK_SET);
  while(1){c=fgetc(fp); if(c!='\n'&&c!=EOF) ct++; else break;}
  fseek(fp,0,SEEK_SET);
  return ct;
}

size_t wordcount(FILE *fp, size_t *linecapa){
  char *line=malloc(*linecapa*sizeof(char)); double num; size_t len,ct=0;
  fseek(fp,0,SEEK_SET);
  getline(&line,linecapa,fp); len=strlen(line);
  fseek(fp,0,SEEK_SET);
  while(fscanf(fp,"%lf",&num)){if(ftell(fp)>=len){break;}ct++;}
  free(line);
  return ct;
}

size_t linecount(FILE *fp, size_t *linecapa){
  char *line=malloc(*linecapa*sizeof(char)); size_t ct=0;
  fseek(fp,0,SEEK_SET);
  while(getline(&line,linecapa,fp)>=0){ct++;}
  fseek(fp,0,SEEK_SET);
  free(line);
  return ct;
}

double *read2dcm(int *nr, int *nc, const char *filename){
  double *a; int n,r,c,i; int sd=sizeof(double); size_t capa;
  FILE *fp=fopen(filename,"r"); if(!fp) goto err01;
  capa=charcount(fp);
  c=wordcount(fp,&capa);
  r=linecount(fp,&capa); n=r*c;
  a=(double*) malloc(n*sd);

  switch(c){
    case 2:  for(i=0;i<r;i++) fscanf(fp,"%lf%lf",   a+c*i,a+c*i+1);          break;
    case 3:  for(i=0;i<r;i++) fscanf(fp,"%lf%lf%lf",a+c*i,a+c*i+1,a+c*i+2);  break;
    default: for(i=0;i<n;i++) fscanf(fp,"%lf",a+i);
  } *nr=r; *nc=c; fclose(fp);

  return a;
  err01: fprintf(stderr,"\n\n  File '%s': Not Found.\n\n",filename); exit(EXIT_FAILURE);
}

