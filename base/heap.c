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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

#define LAST (*heap)
#define ROOT 1

static void swap(int *a, int *b){ int tmp=*a; *a=*b; *b=tmp; }

static void shiftup(int *heap, const double *keys, int i){ 
  int p; assert(i>=1&&i<=LAST); 
  for(;(p=floor(i/2))&&(keys[heap[p]]>keys[heap[i]]);i=p){ swap(heap+p,heap+i); }
}

static void shiftdn(int *heap, const double *keys, int i){
  int c,cl,cr; assert(i);

  for(;(cl=2*i)<=LAST;i=c){ cr=cl+1;
    c=(LAST==cl)?cl:(keys[heap[cl]]<keys[heap[cr]]?cl:cr);
    if(keys[heap[i]]>keys[heap[c]]) swap(heap+c,heap+i);
    else break;
  }
}

static void heap_find(int *i, const int *heap, int node){
  int j;
  for(j=1;j<=LAST&&heap[j]!=node;j++){}
  *i=j>LAST?0:j;
}

void heap_init(int *heap){LAST=0;}

void heap_insert(int *heap, double *keys, int node, double key){
  heap[++LAST]=node;
  keys[node]=key;
  shiftup(heap,keys,LAST);
}

void heap_extract(int *node, int *heap, const double *keys){
  *node=heap[ROOT];
  heap[ROOT]=heap[LAST--];
  shiftdn(heap,keys,ROOT);
}

void heap_downkey(int *heap, double *keys, int node, double key){
  int i; assert(key<keys[node]);
  keys[node]=key;
  heap_find(&i,heap,node); assert(i); //  heap must contain 'node' (=heap[i]).
  shiftup(heap,keys,i);
}  

void print_heap(const int *heap, const double *keys, int mode){
  int i,j,d=1; int L,R;
  for(i=1;i<=LAST;i++){
    if(i%d==0) { d++;
      L=pow(2,d-2);
      R=pow(2,d-1)-1;
      R=LAST<R?LAST:R;
      if(R<=LAST){ 
        for(j=L;j<=R;j++){
          if(mode) printf("%d%c",   heap[j],      (j==R)?'\n':' ');
          else     printf("%.1lf%c",keys[heap[j]],(j==R)?'\n':' ');
        }
      } 
    }
  } printf("\n");

}

