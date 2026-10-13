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
#include<errno.h>
#include<limits.h>
#include<stdint.h>
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

static int separator(int c){return c==','||c==' '||c=='\t'||c=='\r'||c=='\n';}

static int readline2d(FILE *fp, char **line, size_t *capa, size_t *len){
  size_t next; char *tmp;
  if(!*line){*capa=4096;*line=malloc(*capa);if(!*line)return -1;}
  *len=0;**line='\0';
  while(1){
    if(*capa-*len<2){
      if(*capa>SIZE_MAX/2)return -1; next=2*(*capa); if(next>INT_MAX)return -1;
      tmp=realloc(*line,next);if(!tmp)return -1;*line=tmp;*capa=next;
    }
    if(!fgets(*line+*len,(int)(*capa-*len),fp))return *len?1:0;
    *len+=strlen(*line+*len);
    if(*len&&(*line)[*len-1]=='\n')return 1;
    if(feof(fp))return 1;
  }
}

double *read2dcm(int *nr, int *nc, const char *filename){
  FILE *fp; char *line=NULL,*p,*end; double *a=NULL,*tmp,val;
  size_t lcapa=0,len=0,lnum=0,row=0,col=0,rcol,num=0,capa=0,next; int status,invalid_suffix,matrix_too_large;
  *nr=*nc=0; fp=fopen(filename,"r"); if(!fp)goto err01;

  while((status=readline2d(fp,&line,&lcapa,&len))>0){
    lnum++;p=line;while(separator((unsigned char)*p))p++;
    if(*p=='\0'||*p=='#') continue;
    for(rcol=0;;rcol++){
      while(separator((unsigned char)*p))p++;
      if(*p=='\0') break;
      errno=0;val=strtod(p,&end); invalid_suffix=*end&&!separator((unsigned char)*end);
      if(end==p)        {fprintf(stderr,"\n\n  File '%s', line %zu: Invalid numeric value.\n\n",            filename,lnum);goto err02;}
      if(errno==ERANGE) {fprintf(stderr,"\n\n  File '%s', line %zu: Numeric value is out of range.\n\n",    filename,lnum);goto err02;}
      if(!isfinite(val)){fprintf(stderr,"\n\n  File '%s', line %zu: NaN and Inf are not permitted.\n\n",    filename,lnum);goto err02;}
      if(invalid_suffix){fprintf(stderr,"\n\n  File '%s', line %zu: Invalid characters after a number.\n\n",filename,lnum);goto err02;
      }
      if(num==capa){
        if(capa>SIZE_MAX/2) goto err03; next=capa?2*capa:65536;
        if(next>SIZE_MAX/sizeof(double)) goto err03;
        tmp=realloc(a,next*sizeof(double));if(!tmp) goto err03;a=tmp;capa=next;
      }
      a[num++]=val;p=end;
    }
    if(!rcol) continue;
    if(!row) col=rcol;
    else if(rcol!=col){
      fprintf(stderr,"\n\n  File '%s', line %zu: Expected %zu values but found %zu.\n\n",filename,lnum,col,rcol);goto err02;
    }
    row++;
  }
  if(status<0) goto err03; matrix_too_large=row>INT_MAX||col>INT_MAX;
  if(ferror(fp))      {fprintf(stderr,"\n\n  Failed while reading file '%s'.\n\n",  filename); goto err02;}
  if(!row)            {fprintf(stderr,"\n\n  File '%s': No numeric data.\n\n",      filename); goto err02;}
  if(matrix_too_large){fprintf(stderr,"\n\n  Matrix in file '%s' is too large.\n\n",filename); goto err02;}
  tmp=realloc(a,num*sizeof(double));if(tmp)a=tmp;
  free(line);fclose(fp);*nr=(int)row;*nc=(int)col;return a;

  err03: fprintf(stderr,"\n\n  read2dcm: Memory allocation failed while reading file '%s'.\n\n",filename);
  err02: free(a);free(line);fclose(fp);exit(EXIT_FAILURE);
  err01: fprintf(stderr,"\n\n  File '%s': Not Found.\n\n",filename);exit(EXIT_FAILURE);
}

void save_variable(const char *prefix, const char *suffix,const double *var, int D, int J, char *fmt, int trans){
  int d,j; char fn[256],str[256]; FILE *fp;
  strcpy(fn,prefix); strcpy(str,fmt);
  strcat(fn,suffix); strcat(str,"%c");
  if(!(fp=fopen(fn,"w"))){fprintf(stderr,"\n  Error: Could not open file '%s' for writing.\n\n",fn); return;}
  if( trans)for(j=0;j<J;j++)for(d=0;d<D;d++) fprintf(fp,str,var[d+D*j],d==D-1?'\n':'\t');
  if(!trans)for(d=0;d<D;d++)for(j=0;j<J;j++) fprintf(fp,str,var[d+D*j],j==J-1?'\n':'\t');
  fclose(fp);
}
