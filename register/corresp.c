// Copyright (c) 2018-2020 Osamu Hirose
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
#include<errno.h>
#include<limits.h>
#include<stdint.h>
#include<string.h>
#include<math.h>
#ifdef USE_OPENMP
#include<omp.h>
#endif
#include"../base/misc.h"
#include"../base/stat.h"
#include"../base/kdtree.h"
#include"bcpd.h"
#include"corresp.h"

#define SQ(x) ((x)*(x))

/* Return allocated functional bandwidths hf[Df]; NULL for ordinary BCPD. */
double *compute_hf(const double *e, pwsz sz, pwpm pm){
  int d,Df=sz.Df; double zeta=Df?pm.eta*(sz.D/(double)Df):0; double *hf=NULL;

  if(sz.D<=0||Df<0||!(zeta>=0)||!isfinite(zeta)) goto err01;
  if(!Df) return NULL;
  if(!e) goto err01;
  hf=calloc(Df,sizeof(double)); if(!hf) goto err02;
  for(d=0;d<Df;d++){
    if(zeta>0&&(!(e[d]>0)||!isfinite(e[d]))) goto err01;
    hf[d]=zeta>0?e[d]/sqrt(zeta):INFINITY;
    if(zeta>0&&(!(hf[d]>0)||!isfinite(hf[d]))) goto err01;
  }
  return hf;
  err01: free(hf);fprintf(stderr,"ERROR: Invalid residual scales or feature weight in correspondence export.\n"); exit(EXIT_FAILURE);
  err02: free(hf);fprintf(stderr,"ERROR: Out of memory preparing correspondence bandwidths.\n"); exit(EXIT_FAILURE);
}

/* Return allocated log outlier contributions lc[N0]; */
double *compute_lc(const double *X, const double *fx, const double *fx0, const double *e, double r, int N0, pwsz sz, pwpm pm){
  int d,n,D=sz.D,Df=sz.Df,N=sz.N; double zeta=Df?pm.eta*(D/(double)Df):0;
  double val; double *lc=NULL,*mf=NULL,*vf=NULL;

  if(D<=0||Df<0||N<=0||N0<=0||!X||!(r>0)||!isfinite(r)||!(zeta>=0)||!isfinite(zeta)) goto err01;
  if(Df&&(!fx||!fx0||!e)) goto err01;
  if(zeta>0) for(d=0;d<Df;d++) if(!(e[d]>0)||!isfinite(e[d])) goto err01;
  if(!(pm.omg>=0&&pm.omg<1)) goto err01;
  lc=calloc(N0,sizeof(double)); if(!lc) goto err02;
  if(pm.omg==0){for(n=0;n<N0;n++) lc[n]=-INFINITY;return lc;}
  val=log(pm.omg)-log1p(-pm.omg)-lnvolume(X,D,N)-lnidet(r,D);
  if(zeta>0){
    val-=zeta*lnddet(e,Df);
    if(Df<=10) val-=zeta*lnvolume(fx,Df,N);
    else{
      mf=calloc(Df,sizeof(double)); vf=calloc(Df,sizeof(double));
      if(!mf||!vf) goto err02;
      mvmean(mf,fx,Df,N); mvsdev(vf,fx,Df,N);
    }
  }
  if(!isfinite(val)) goto err01;
  for(n=0;n<N0;n++){
    lc[n]=val+((zeta>0&&Df>10)?zeta*lnnormd(fx0+Df*n,mf,vf,Df):0);
    if(isnan(lc[n])||lc[n]==INFINITY) goto err01;
  }
  free(mf);free(vf);return lc;
  err01: free(lc);free(mf);free(vf);fprintf(stderr,"ERROR: Invalid residual scales or outlier model in correspondence export.\n"); exit(EXIT_FAILURE);
  err02: free(lc);free(mf);free(vf);fprintf(stderr,"ERROR: Out of memory preparing correspondence outliers.\n"); exit(EXIT_FAILURE);
}

/* Extend learned source masses by nearest-sample cells in the original, normalized source frame. */
double *extend_a(
    const double *Y0, const double *Y, const double *a, const int *U,
    int D, int M0, int M, int learned
  ){
  int i,m; int *T=NULL,*rep=NULL,*id=NULL,*count=NULL; double *ar=NULL,*a0,dist;
  int si=sizeof(int),sd=sizeof(double);

  a0=calloc(M0,sd); if(!a0) goto err01;
  if(!learned){for(m=0;m<M0;m++) a0[m]=1.0/M0;return a0;}
  T=kdtree_build(Y,D,M);
  rep=calloc(M,si); ar=calloc(M,sd);
  id=calloc(M0,si); count=calloc(M,si);
  if(!T||!rep||!ar||!id||!count) goto err01;
  /* rep maps samples to representatives; ar holds their aggregated masses. */
  for(m=0;m<M0;m++) id[m]=-1;
  for(i=0;i<M;i++){
    if(U[i]<0||U[i]>=M0||!(a[i]>=0)||!isfinite(a[i])) goto err02;
    if(id[U[i]]<0) id[U[i]]=i;
    rep[i]=id[U[i]];ar[rep[i]]+=a[i];
  }
  for(m=0;m<M0;m++){
    if(id[m]<0){nnsearch(&i,&dist,Y0+D*m,Y,T,D,M);if(i<0) goto err02;id[m]=rep[i];}
    count[id[m]]++;
  }
  for(m=0;m<M0;m++){i=id[m];a0[m]=ar[i]/count[i];}
  free(T);free(rep);free(ar);free(id);free(count);return a0;
  err01: fprintf(stderr,"ERROR: Out of memory extending correspondence weights.\n"); exit(EXIT_FAILURE);
  err02: fprintf(stderr,"ERROR: Invalid source sampling map in correspondence export.\n"); exit(EXIT_FAILURE);
}

static FILE *corresp_open(const char *prefix, const char *suffix, const char *header){
  char fn[512]; FILE *fp; int n=snprintf(fn,sizeof(fn),"%s%s",prefix,suffix);
  if(n<0||(size_t)n>=sizeof(fn)){fprintf(stderr,"ERROR: Correspondence filename is too long.\n");exit(EXIT_FAILURE);}
  fp=fopen(fn,"w");
  if(!fp){fprintf(stderr,"ERROR: Cannot open correspondence file '%s': %s\n",fn,strerror(errno));exit(EXIT_FAILURE);}
  if(fputs(header,fp)==EOF){fprintf(stderr,"ERROR: Cannot write correspondence file '%s'.\n",fn);exit(EXIT_FAILURE);}
  return fp;
}

static void corresp_close(FILE *fp){
  int err; if(!fp) return;
  err=ferror(fp); if(fclose(fp)) err=1;
  if(err){fprintf(stderr,"ERROR: Cannot finish writing correspondence output.\n");exit(EXIT_FAILURE);}
}

static void corresp_error(int err, int n){
  if(err==1) fprintf(stderr,"ERROR: Invalid input to correspondence export.\n");
  if(err==2) fprintf(stderr,"ERROR: Out of memory in correspondence export.\n");
  if(err==3) fprintf(stderr,"ERROR: No finite correspondence score for target %d.\n",n+1);
  exit(EXIT_FAILURE);
}

/* Approximate correspondences at the returned fitted deformation. */
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
  ){
  int d,m,D,Df,M,N,nt=1; int *T; double rad,*la,*x=NULL,*wc=NULL;
  int si=sizeof(int),sd=sizeof(double); FILE *fpP=NULL,*fpc=NULL,*fpe=NULL; int nzero=0;

  D=sz.D; Df=sz.Df; M=sz.M; N=sz.N; rad=pm.dlt*r;
  if(D<=0||M<=0||N<=0||M>INT_MAX/6||!(r>0)||!isfinite(r)||!(rad>0)||isnan(rad)) goto err01;
  if(Df&&(!fx||!fy||!hf)) goto err01;
  if(pm.opt&PW_OPT_SAVEP) fpP=corresp_open(prefix,"P.txt","[n]\t[m]\t[probability]\n");
  if(pm.opt&PW_OPT_SAVEE) fpe=corresp_open(prefix,"e.txt","[n]\t[m]\t[probability]\n");
  if(pm.opt&PW_OPT_SAVEC) fpc=corresp_open(prefix,"c.txt","[n]\t[1/0]\n");

  T=kdtree_build(y,D,M);
  la=calloc(M,sd);
  if(!T||!la) goto err02;
  if((size_t)M>SIZE_MAX/(size_t)D/sd) goto err02;
  if(pm.opt&PW_OPT_SAVEX){
    x=calloc((size_t)D*M,sd); wc=calloc(M,sd);
    if(!x||!wc) goto err02;
  }
  /* Ordered text exports are faster serially. Bound parallel scratch to 512 MiB. */
  #ifdef USE_OPENMP
  if(x&&!fpP&&!fpc&&!fpe){
    nt=(int)(536870912. / ((double)M*((D+2.)*sd+si)));
    if(nt<1) nt=1;
    if(nt>omp_get_max_threads()) nt=omp_get_max_threads();
    if(nt>N) nt=N;
  }
  #endif
  /* Cache log mixture coefficients independently of the tree workspace. */
  for(m=0;m<M;m++){if(!(a[m]>=0)||!isfinite(a[m])) goto err01;la[m]=a[m]>0?log(a[m]):-INFINITY;}
  #pragma omp parallel num_threads(nt) private(d,m)
  {
    int i,n,top,ct,mmax,S[MAXTREEDEPTH]; int *l=calloc(M,si);
    double c,val,max,min,*p=calloc(M,sd),*xs=x,*ws=wc;
    if(x&&nt>1){xs=calloc((size_t)D*M,sd);ws=calloc(M,sd);}
    if(!l||!p||(x&&(!xs||!ws))) corresp_error(2,0);
    #pragma omp for schedule(dynamic,64)
    for(n=0;n<N;n++){
      /* Radius search, with nearest-neighbor fallback. */
      top=ct=0;do{eballsearch_next(&m,S,&top,X+D*n,rad,y,T,D,M);if(m>=0)l[ct++]=m;}while(top);
      if(!ct){nnsearch(&m,&min,X+D*n,y,T,D,M);if(m<0) corresp_error(1,n);l[ct++]=m;}
      /* Normalize log scores to avoid kernel underflow. */
      max=lc[n];mmax=0;
      if(isnan(max)||max==INFINITY) corresp_error(1,n);
      for(i=0;i<ct;i++){m=l[i];val=0;
        for(d=0;d<D; d++) val+=SQ((y[d+D*m]-X[d+D*n])/r);
        for(d=0;d<Df;d++) val+=SQ((fy[d+Df*m]-fx[d+Df*n])/hf[d]);
        p[i]=la[m]-0.5*val;
        if(isnan(p[i])||p[i]==INFINITY) corresp_error(1,n);
        if(p[i]>max) max=p[i];
        if(p[i]>p[mmax]) mmax=i;
      }
      if(!isfinite(max)) corresp_error(3,n);
      c=exp(lc[n]-max);val=c;
      for(i=0;i<ct;i++){p[i]=exp(p[i]-max);val+=p[i];}
      c/=val;for(i=0;i<ct;i++) p[i]/=val;
      /* Accumulate all probabilities, including those below the export threshold. */
      if(xs) for(i=0;i<ct;i++){
        m=l[i];ws[m]+=p[i];
        for(d=0;d<D;d++) xs[d+D*m]+=p[i]*X[d+D*n];
      }
      /* e reports the best source even when c marks it as an outlier. */
      if(fpP){for(i=0;i<ct;i++)if(p[i]>1.0/M){m=l[i];fprintf(fpP,"%d\t%d\t%.17g\n",n+1,m+1,p[i]);}}
      if(fpe){fprintf(fpe,"%d\t%d\t%.17g\n",n+1,l[mmax]+1,p[mmax]);}
      if(fpc){fprintf(fpc,"%d\t%d\n",n+1,p[mmax]>c?1:0);}
    }
    /* Merge each worker once; search and probability updates need no locks. */
    if(x&&nt>1){
      #pragma omp critical
      {for(m=0;m<M;m++){ wc[m]+=ws[m]; for(d=0;d<D;d++) x[d+D*m]+=xs[d+D*m];}}
      free(xs);free(ws);
    }
    free(l);free(p);
  }
  if(x){
    for(m=0;m<M;m++){
      if(wc[m]>0){for(d=0;d<D;d++) x[d+D*m]/=wc[m];}
      else{for(d=0;d<D;d++) x[d+D*m]=y[d+D*m];nzero++;}
    }
    if(nzero&&!(pm.opt&PW_OPT_QUIET))
      fprintf(stderr,"  Correspondence export: %d zero-mass source points use y for x.\n",nzero);
  }
  corresp_close(fpP);corresp_close(fpc);corresp_close(fpe);
  free(wc);free(la);free(T);return x;
  err01: corresp_error(1,0);return NULL;
  err02: corresp_error(2,0);return NULL;
}
