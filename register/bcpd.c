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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>
#include<time.h>
#include<string.h>
#include<sys/time.h>
#include"../base/misc.h"
#include"../base/lapack.h"
#include"../base/kdtree.h"
#include"../base/kernel.h"
#include"../base/gaussprod.h"
#include"../base/gramdecomp.h"
#include"bcpd.h"

int dgemm_ (char *transA,char *transB,int *M,int*N,int*K,double*alpha,double*A,int*lda,double*B,int*ldb,double*beta,double*C,int*ldc);

#define SQ(x) ((x)*(x))
double digamma(double x);
double (*kernel[4])(const double*,const double*,int,double)={gauss,imquad,rational,laplace};
char   spinner[16][16]={
  "-    \0", "\\    \0",  "|    \0", "/    \0",
  "-    \0", " -   \0",   "  -  \0", "   - \0",
  "    -\0", "    \\\0",  "    |\0", "    /\0",
  "    -\0", "   - \0",   "  -  \0", " -   \0",
};

static void print_status(int lp, double Np, double sigma, double diff, double conv, int opt, int J, int local){
  char *acc=(local&&(opt&PW_OPT_LOCAL))?"KDtree ":(J?"Nystrom":"None   ");
  if(opt&PW_OPT_QUIET) return;
  if(opt&PW_OPT_HISTO){
    if(lp) fprintf(stderr,"  loop=%.3d  acc.P=[%s]  N_hat=%lf  sigma=%lf  diff=%lf\n",lp+1,acc,Np,sigma,diff);
    else   fprintf(stderr,"  loop=%.3d  acc.P=[%s]  N_hat=%lf  sigma=%lf          \n",lp+1,acc,Np,sigma);
  }
  else {
    if(lp) fprintf(stderr,"\033[8A");
    fprintf(stderr,"  VB optimization: -|%s|-\n", spinner[lp%16]);
    fprintf(stderr,"    loop   =  %d\n",   lp+1);
    fprintf(stderr,"    Np     =  %lf\n",  Np);
    fprintf(stderr,"    sigma  =  %lf\n",  sigma);
    fprintf(stderr,"    diff   =  "); if(lp) fprintf(stderr,"%lf\n",diff); else fprintf(stderr,"\n");
    fprintf(stderr,"    conv   =  %lf\n",  conv);
    fprintf(stderr,"    acc.P  =  %s\n\n", acc);
  }
}

static void lowrank_encode(
  double         *B,
  double         *wd, /* size: MD+2KK+(D)MK (DMK if ssm) */
  const double   *E,
  const double   *w,
  const double   *s,
  const double   *r,
  const double   *L,
  const double   *Q,
  int             D,
  int             M,
  int             K,
  double          lmd,
  int             ssm
  ){

  int d,i,j,k,m; double *A,*C,*S,*W; char uplo='U'; int D0=D;
  int info; double val,cc; int sd=0;

  if(ssm){M=D*M;D=1;} /* ssm order: xyzxyz */

  /* allocation */
  W=wd+sd;sd+=D*M; C=wd+sd;sd+=M*K;
  A=wd+sd;sd+=K*K; S=wd+sd;sd+=K*K;

  cc=lmd*SQ(*r/(*s));
  for(m=0;m<M;m++)for(k=0;k<K;k++) C[m+M*k]=w[ssm?(m/D0):m]*Q[m+M*k];
  for(k=0;k<K;k++)for(d=0;d<D;d++){B[k+K*d]=0;for(m=0;m<M;m++) B[k+K*d]+=C[m+M*k]*E[d+D*m];}
  #pragma omp parallel for private (j) private (m) private(val)
  for(i=0;i<K;i++)for(j=0;j<K;j++){val=0;for(m=0;m<M;m++){val+=Q[m+M*i]*C[m+M*j];} A[i+K*j]=val;}
  for(i=0;i<K;i++)for(j=0;j<K;j++) S[i+K*j]=A[i+K*j];
  for(k=0;k<K;k++) A[k+K*k]+=cc/L[k];
  dpotrf_(&uplo,&K,A,&K,&info);         if(info!=0){goto err02;}
  dpotrs_(&uplo,&K,&D,A,&K,B,&K,&info); if(info!=0){goto err03;}
  dpotrs_(&uplo,&K,&K,A,&K,S,&K,&info); if(info!=0){goto err04;}
  for(i=0;i<K;i++)for(j=0;j<K;j++) A[i+K*j]=L[i]*((i==j?1:0)-S[j+K*i]);
  for(m=0;m<M;m++)for(d=0;d<D;d++) W[m+M*d]=w[ssm?(m/D0):m]*E[d+D*m];
  for(m=0;m<M;m++)for(d=0;d<D;d++)for(k=0;k<K;k++) W[m+M*d]-=C[m+M*k]*B[k+K*d];
  for(m=0;m<M;m++)for(d=0;d<D;d++) W[m+M*d]/=cc;
  for(k=0;k<K;k++)for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=Q[m+M*k]*W[m+M*d];} B[k+K*d]=val*L[k];}

  return;
  err02: printf("ERROR: The Cholesky factorization of A failed at the update of v. Retry.\n");     exit(EXIT_FAILURE);
  err03: printf("ERROR: Solving linear equations AX=B failed at the update of v.\n");              exit(EXIT_FAILURE);
  err04: printf("ERROR: Solving linear equations AX=S failed at the update of v.\n");              exit(EXIT_FAILURE);
}

static void lowrank_decode(
  double         *u,
  double         *v,
  const double   *Y,
  const double   *Q,
  const double   *B,
  int             D,
  int             M,
  int             K,
  int             ssm
  ){
  int d,k,m; double val;
  if(ssm){M=D*M;D=1;}
  for(d=0;d<D;d++)for(m=0;m<M;m++){val=0;for(k=0;k<K;k++){val+=Q[m+M*k]*B[k+K*d];} v[d+D*m]=u[d+D*m]=val; u[d+D*m]+=Y[d+D*m];}
}

static void invlinear(double *z, const double *y, int D, int M, const double *s, const double *R, const double *t, char type, double *w, int *ipiv){
  int d,i,j,m; int info; assert(type=='s'||type=='a'); /* similarity or affine */

  if(type=='s')for(m=0;m<M;m++)for(d=0;d<D;d++){z[d+D*m]=0;for(i=0;i<D;i++){z[d+D*m]+=R[i+D*d]*(y[i+D*m]-t[i]);} z[d+D*m]/=(*s);}
  if(type=='a'){ assert(*s==1.0);
    for(i=0;i<D;i++)for(j=0;j<D;j++) w[i+D*j]=R[i+D*j];
    for(m=0;m<M;m++)for(d=0;d<D;d++) z[d+D*m]=y[d+D*m]-t[d];
    dgesv_(&D,&M,w,&D,ipiv,z,&D,&info);assert(!info);
  }
}

int bcpd(
  double        *   x,    /*  O  | DM x 1 (+nlp) | aligned target shape     */
  double        *   y,    /*  O  | DM x 1 (+nlp) | deformed source shape    */
  double        *   u,    /*  O  | DM x 1        | normalized def. shape    */
  double        *   v,    /*  O  | DM x 1        | displacement vectors     */
  double        *   w,    /*  O  |  M x 1        | #matches for each m      */
  double        *   a,    /*  O  |  M x 1        | mixing coefficients      */
  double        *   sgm,  /*  O  |  M x 1        | posterior covariance     */
  double        *   s,    /*  O  |    1          | scale factor             */
  double        *   R,    /*  O  |  D x D        | rotation matrix          */
  double        *   t,    /*  O  |    D          | translation vector       */
  double        *   r,    /*  O  |    1          | residual s.d.            */
  double        *   e,    /*  O  |    Df         | residual s.d. of fx, fy  */
  double        *   Np,   /*  O  |    1          | #matched points (est'd)  */
  double        *   pf,   /*  O  | nlp x 3       | comp. time (r/c) & sigma */
  double        *   wd,   /*  W  |    *          | working memory (double)  */
  int           *   wi,   /*  W  |    *          | working memory (int)     */
  const double  *   X,    /*  I  | DN x 1        | target point set         */
  const double  *   Y,    /*  I  | DM x 1        | source point set         */
  const double  *   fx,   /*  I  |  N x 1        | function value of x      */
  const double  *   fy,   /*  I  |  M x 1        | function value of y      */
  const double  *   LQ,   /*  I  | K + M x K     | only for geodesic kernel */
  const pwsz        sz,   /*  I  |               | D, Df, M, N, K, J,       */
  const pwpm        pm    /*  I  |               | tuning parameters        */
  ){

  double c,cc,val,val1,val2,c1,c2; double vol,diff,rold=1e100,reg=1e-20,dlt,lim; int flg,local; int mtd;
  double *B/*KD*/,*E/*MD*/,*Q/*MK*/,*L/*K*/,*G/*MM*/,*G1/*MM*/,*G2/*MM*/;
  double *b/*M*/,*q/*N*/,*f/*N*/,*PX/*DM*/,*sx/*D*M*nlp*/,*sy/*D*M*nlp*/,*ix/*DM*/,*lw/*DD*/;
  double *xb/*D*/,*ub/*D*/,*phi/*DD*/,*psi/*DD*/,*Sxu/*DD*/,*dS/*D*/,*wk/*10D*/; struct timeval tick;
  double *wld,*wdd/*K(K+11)*/,*wgd/*M+N+J+JJ*/; int *Tx,*Ty/*1+3N*/,*wdi/*M*/,*wgi/*M+N*/; double bet;
  int d,i,j,m,n,lp,D,Df,K,M,N,J; int *ipiv,info,lwork1; char jobz='V',uplo='U'; double tr=0; int max;
  double omg,lmd,kpa,cnv; int nlp,opt; int sd=0,si=0; int T=pm.opt&PW_OPT_LOCAL,db=pm.opt&PW_OPT_DBIAS;
  double clc,cps=CLOCKS_PER_SEC,*treal=pf,*tcpu=pf+pm.nlp,*rprog=pf+2*pm.nlp; int pflog=pm.opt&PW_OPT_PFLOG;
  double *Pfx; int ssm=strlen(pm.fn[COV_LQ]); int type=(pm.opt&PW_OPT_AFFIN)?'a':'s';

  /* record base time */
  if(pflog){ gettimeofday(&tick,NULL); clc=(double)clock()/cps;
    for(lp=0;lp<pm.nlp;lp++) treal[lp]=-(tick.tv_sec+tick.tv_usec/1e6);
    for(lp=0;lp<pm.nlp;lp++) tcpu [lp]=-clc;
  }
  /* alias: parameter & size */
  omg=pm.omg; kpa=pm.kpa; cnv=pm.cnv; D=sz.D; K=sz.K; mtd=MAXTREEDEPTH;
  lmd=pm.lmd; dlt=pm.dlt; nlp=pm.nlp; M=sz.M; J=sz.J; lwork1=10*D;
  bet=pm.bet; lim=pm.lim; opt=pm.opt; N=sz.N; max=M+N;/*max=M>N?M:N*/; Df=sz.Df;
  /*---------------------------------------------------------------o
  |   alias: memory                                                |
  o---------------------------------------------------------------*/
  /* common: double-> 4M+2N+D(5M+N+13D+3), int->D                 */
  b=wd+sd; sd+=M; lw=wd+sd; sd+=D*D; phi=wd+sd; sd+=D*D; dS=wd+sd; sd+=D;
  /*-----------*/ PX=wd+sd; sd+=M*D; psi=wd+sd; sd+=D*D; xb=wd+sd; sd+=D;
  f=wd+sd; sd+=N; /*--------------*/ Sxu=wd+sd; sd+=D*D; ub=wd+sd; sd+=D;
  q=wd+sd; sd+=N; E =wd+sd; sd+=D*M; ix =wd+sd; sd+=D*M; wk=wd+sd; sd+=lwork1;
  ipiv=wi+si; si+=D; Pfx=Df?wd+sd:NULL; sd+=Df?Df*M:0;

  if(pm.opt&PW_OPT_NONRG){G=G1=G2=NULL; goto skip_G_alloc;}
  /* rank restriction: double-> K x (2M+3K+D+12), int->M */
  /*------------------------*/ G  =K?NULL:wd+sd; sd+=K?0:M*M;
  B=K?wd+sd:NULL; sd+=K?K*D:0; G1 =K?NULL:wd+sd; sd+=K?0:M*M;
  /*------------------------*/ G2 =K?NULL:wd+sd; sd+=K?0:M*M;
  if(!LQ){
    Q=K?wd+sd:NULL; sd+=M*K*(ssm?D:1);
    L=K?wd+sd:NULL; sd+=K?K:0;
  } else  L=Q=NULL;
  wdd=K?wd+sd:NULL; sd+=K?K*(K+11):0;
  wdi=K?wi+si:NULL; si+=K?M:0;
  wld=K?wd+sd:NULL; sd+=M*K*(ssm?D:1)+(K?M*D+2*K*K:0);
  skip_G_alloc:
  /* gaussprod */
  wgd=(J||T)?wd+sd:NULL;
  wgi=(J||T)?wi+si:NULL;
  if(J){sd+=2*J*(1+D+Df+J); si+=M+N;}
  if(T){sd+=2*max;si+=max*(8+mtd);}
  /* tree */
  Tx=T?wi+si:NULL; si+=T?3*max+1:0;
  Ty=T?wi+si:NULL; si+=T?3*max+1:0;
  /* trajectory */
  sx=(pm.opt&PW_OPT_SAVEX)?(x+D*M):NULL;
  sy=(pm.opt&PW_OPT_SAVEY)?(y+D*M):NULL;

  /*---------------------------------------------------------------o
  |   initialization                                               |
  o---------------------------------------------------------------*/
  /* a */
  for(m=0;m<M;m++) a[m]=1/(double)M;
  for(m=0;m<M;m++) sgm[m]=0;
  /* s,R,t */
  *s=1.0;
  for(d=0;d<D;d++)for(i=0;i<D;i++) R[d+D*i]=d==i?1:0;
  for(d=0;d<D;d++) t[d]=0;
  /* y, b */
  for(d=0;d<D;d++)for(m=0;m<M;m++) y[d+D*m]=Y[d+D*m];
  for(m=0;m<M;m++) b[m]=1.0;
  /* r */
  *r=0;
  for(d=0;d<D;d++){ val1=val2=0;
    for(m=0;m<M;m++) *r+=N*SQ(Y[d+D*m]);
    for(n=0;n<N;n++) *r+=M*SQ(X[d+D*n]);
    for(m=0;m<M;m++) val1+=Y[d+D*m];
    for(n=0;n<N;n++) val2+=X[d+D*n];
    *r-=2*val1*val2;
  } *r/=(size_t)M*N*D; *r=sqrt(*r)*pm.gma;
  /* e */
  for(d=0;d<Df;d++){ e[d]=val1=val2=0;
    for(m=0;m<M;m++) e[d]+=N*SQ(fy[d+Df*m]);
    for(n=0;n<N;n++) e[d]+=M*SQ(fx[d+Df*n]);
    for(m=0;m<M;m++) val1+=fy[d+Df*m];
    for(n=0;n<N;n++) val2+=fx[d+Df*n];
    e[d]-=2*val1*val2;
    e[d]/=(size_t)M*N; e[d]=sqrt(e[d])*pm.gma;
  }
  if(pm.opt&PW_OPT_ISOVF){val=0;
    for(d=0;d<Df;d++) val+=SQ(e[d]);
    for(d=0;d<Df;d++) e[d]=sqrt(val/Df);
  }
  /* volume */
  vol=volume(X,D,N)*(Df?volume(fx,Df,N):1.0f);
  /* tree */
  if(T) kdtree(Tx,wgi,wgd,X,D,N);
  /* G */
  if(pm.opt&PW_OPT_NONRG){for(d=0;d<D;d++)for(m=0;m<M;m++){u[d+D*m]=Y[d+D*m];} goto skip_G;}
  if(!K)
    #pragma omp parallel for private (j)
    for(i=0;i<M;i++)for(j=i;j<M;j++) G[i+M*j]=G[j+M*i]=kernel[pm.G](Y+D*i,Y+D*j,D,bet);
  else if(LQ){L=(double*)LQ;Q=(double*)LQ+K;}
  else gramdecomp(Q,L,wdd,wdi,Y,D,M,K,bet,kernel[pm.G]);
  skip_G:
  /* save comp. profile */
  if(pflog){gettimeofday(&tick,NULL);treal[0]+=tick.tv_sec+tick.tv_usec/1e6;tcpu[0]+=clock()/cps;rprog[0]=*r;}

  /*---------------------------------------------------------------o
  |   VB optimization                                              |
  o---------------------------------------------------------------*/
  for(lp=0;lp<nlp;lp++){
    /*---------------------------------------------------------------o
    |   update: x, w                                                 |
    o---------------------------------------------------------------*/
    local=(*r<pm.btn); flg=(local&&(opt&PW_OPT_LOCAL))?GRAM_FLAG_LOCAL:0;
    c=(pow(2.0*M_PI*SQ(*r),0.5*D)*omg)/(vol*(1-omg)); /* c */
    if(Df) c*=pow(2.0*M_PI,0.5*Df)*product(e,Df);
    for(m=0;m<M;m++) b[m]=a[m]*exp(-(D/2)*sgm[m]*SQ((*s)/(*r)));
    gaussprod(q,NULL,NULL,wgd,wgi,y,X,fy,fx,b,Ty,D,Df,M,N,J,*r,e,dlt,lim,flg|GRAM_FLAG_TRANS|GRAM_FLAG_BUILD);/* Kt1 */
    for(n=0;n<N;n++) q[n]=1.0/(q[n]+c); /* q */
    for(n=0;n<N;n++) f[n]=1.0-(q[n]*c); /* f */
    gaussprod(w,PX,Pfx,wgd,wgi,y,X,fy,fx,q,Tx,D,Df,M,N,J,*r,e,dlt,lim,flg);
    for(m=0;m<M;m++){w[m]*=b[m];w[m]=w[m]<reg?reg:w[m];}
    for(d=0;d<D;d++)for(m=0;m<M;m++) PX[m+M*d]*=b[m];
    if(Df) for(d=0;d<Df;d++)for(m=0;m<M;m++) Pfx[m+M*d]*=b[m];
    *Np=0;for(m=0;m<M;m++) *Np+=w[m]; /* Np */
    for(m=0;m<M;m++)for(d=0;d<D;d++) x[d+D*m]=PX[m+M*d]/w[m]; /* x */
    /*---------------------------------------------------------------o
    |   save trajectory                                              |
    o---------------------------------------------------------------*/
    if(opt&PW_OPT_PATHX) for(m=0;m<M;m++)for(d=0;d<D;d++) sx[d+D*m+D*M*lp]=x[d+D*m];
    if(opt&PW_OPT_PATHY) for(m=0;m<M;m++)for(d=0;d<D;d++) sy[d+D*m+D*M*lp]=y[d+D*m];
    /*---------------------------------------------------------------o
    |   update: u, v                                                 |
    o---------------------------------------------------------------*/
    if(pm.opt&PW_OPT_NONRG) goto skip_uv;
    cc=lmd*SQ(*r/(*s)); invlinear(ix,x,D,M,s,R,t,type,lw,ipiv);
    for(m=0;m<M;m++)for(d=0;d<D;d++) E[K?d+D*m:m+M*d]=ix[d+D*m]-Y[d+D*m];
    if(K){/* CASE: low rank  */
      lowrank_encode(B,wld,E,w,s,r,L,Q,D,M,K,lmd,ssm);
      lowrank_decode(u,v,Y,Q,B,D,M,K,ssm);
    }
    else{ /* CASE: full rank */
      #pragma omp parallel for private (j)
      for(i=0;i<M;i++)for(j=i;j<M;j++) G1[i+M*j]=G1[j+M*i]=G[i+M*j]+(i==j?cc/w[i]:0);
      dposv_(&uplo,&M,&D,G1,&M,E,&M,&info); if(info!=0){goto err06;}
      for(m=0;m<M;m++)for(d=0;d<D;d++) u[d+D*m]=Y[d+D*m];
      #pragma omp parallel for private (d) private (i)
      for(m=0;m<M;m++)for(d=0;d<D;d++){u[d+D*m]=Y[d+D*m];for(i=0;i<M;i++)u[d+D*m]+=G[m+M*i]*E[i+M*d];}
      if(db){tr=0;
        for(i=0;i<M;i++)for(j=i;j<M;j++) G2[i+M*j]=G2[j+M*i]=G[i+M*j];
        dpotrs_(&uplo,&M,&M,G1,&M,G2,&M,&info);if(info!=0){goto err07;}
        for(m=0;m<M;m++){sgm[m]=G2[m+M*m]*SQ(*r/(*s))/w[m]; tr+=sgm[m];} tr*=D;
      }
    }
    skip_uv:
    /*---------------------------------------------------------------o
    |   update: a                                                    |
    o---------------------------------------------------------------*/
    if(kpa>ZERO) for(m=0;m<M;m++) a[m]=exp(digamma(kpa+w[m])-digamma(*Np+kpa*M));
    /*---------------------------------------------------------------o
    |   update: s, R, t                                              |
    o---------------------------------------------------------------*/
    if((pm.opt&PW_OPT_NOSIM)&&(!(pm.opt&PW_OPT_AFFIN))) goto skip_T;
    /* xb, yb */
    for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=u[d+D*m]*w[m];} ub[d]=val/(*Np);}
    for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=x[d+D*m]*w[m];} xb[d]=val/(*Np);}
    /* Sxu & Suu */
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(m=0;m<M;m++){val+=w[m]*(x[i+D*m]-xb[i])*(u[j+D*m]-ub[j]);} Sxu[type=='a'?(j+D*i):(i+D*j)]=val;}
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(m=0;m<M;m++){val+=w[m]*(u[i+D*m]-ub[i])*(u[j+D*m]-ub[j]);} R  [i+D*j]=val;}
    /* A */
    if(type=='a'){dposv_(&uplo,&D,&D,R,&D,Sxu,&D,&info);assert(!info); for(i=0;i<D;i++)for(j=0;j<D;j++){R[i+D*j]=Sxu[j+D*i];} goto skip_sR;}
    /* R */
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(d=0;d<D;d++){val+=Sxu[d+D*i]*Sxu[d+D*j];} psi[i+D*j]=val;}
    dsyev_(&jobz,&uplo,&D,psi,&D,dS,wk,&lwork1,&info);
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(d=0;d<D;d++){val+=Sxu[i+D*d]*psi[d+D*j];} phi[i+D*j]=val/sqrt(dS[j]);}
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(d=0;d<D;d++){val+=phi[i+D*d]*psi[j+D*d];} R[i+D*j]=val;}
    dgetrf_(&D,&D,R,&D,ipiv,&info);if(info!=0){goto err08;}
    val=1.0;for(d=0;d<D;d++){val*=R[d+D*d];if(d+1!=ipiv[d])val*=-1.0;}
    if(val<0)for(d=0;d<D;d++)phi[d+D*(D-1)]*=-1;
    for(i=0;i<D;i++)for(j=0;j<D;j++){val=0;for(d=0;d<D;d++){val+=phi[i+D*d]*psi[j+D*d];} R[i+D*j]=val;}
    /* s */
    c1=c2=0; if(db) c2+=tr/(*Np);
    for(i=0;i<D;i++)for(j=0;j<D;j++) c1+=R[i+D*j]*Sxu[i+D*j];
    for(d=0;d<D;d++)for(m=0;m<M;m++) c2+=w[m]*(SQ(u[d+D*m]-ub[d]));
    *s=(pm.opt&PW_OPT_NOSCL)?1.0:(c1/c2);
    skip_sR:
    /* t */
    for(d=0;d<D;d++){val=0;for(i=0;i<D;i++)val+=R[d+D*i]*ub[i];t[d]=xb[d]-(*s)*val;}
    skip_T:
    /* y */
    for(d=0;d<D;d++)for(m=0;m<M;m++){val=0;for(i=0;i<D;i++){val+=R[d+D*i]*u[i+D*m];} y[d+D*m]=(*s)*val+t[d];}
    /*---------------------------------------------------------------o
    |   update: residual                                             |
    o---------------------------------------------------------------*/
    rold=*r; *r=0; if(db) *r+=SQ(*s)*tr;
    for(n=0;n<N;n++)for(d=0;d<D;d++) *r+=SQ(X[d+D*n])*f[n];
    for(m=0;m<M;m++)for(d=0;d<D;d++) *r+=SQ(y[d+D*m])*w[m];
    for(m=0;m<M;m++)for(d=0;d<D;d++) *r-=2*PX[m+M*d]*y[d+D*m];
    if(*Np<0){goto err05;}
    *r/=(*Np)*D; *r=fabs(*r); *r=sqrt(*r);
    diff=fabs(rold-*r);
    print_status(lp,*Np,*r,diff,cnv,opt,J,local);
    if(lp>pm.llp&&(*r<cnv||diff<cnv)){
      if(kpa>ZERO)for(m=0;m<M;m++) a[m]=(kpa+w[m])/(kpa*M+(*Np));
      break;
    }
    /* update e */
    for(d=0;d<Df;d++){ e[d]=0;
      for(n=0;n<N;n++) e[d]+=SQ(fx[d+Df*n])*f[n];
      for(m=0;m<M;m++) e[d]+=SQ(fy[d+Df*m])*w[m];
      for(m=0;m<M;m++) e[d]-=2*Pfx[m+M*d]*fy[d+Df*m];
      e[d]/=(*Np); e[d]=fabs(e[d]); e[d]=sqrt(e[d]); e[d]+=pm.fpv;
    }
    if(pm.opt&PW_OPT_ISOVF){val=0;
      for(d=0;d<Df;d++) val+=SQ(e[d]);
      for(d=0;d<Df;d++) e[d]=sqrt(val/Df);
    }
    /* save comp. profile */
    if(pflog&&lp+1<nlp){gettimeofday(&tick,NULL);treal[lp+1]+=tick.tv_sec+tick.tv_usec/1e6;tcpu[lp+1]+=clock()/cps;rprog[lp+1]=*r;}
  }

  return lp;

  err05: printf("ERROR: Np became negative. Abort.\n");                                            exit(EXIT_FAILURE);
  err06: printf("ERROR: The Cholesky factorization of G failed at the update of v. Retry.\n");     exit(EXIT_FAILURE);
  err07: printf("ERROR: Inversion of G failed at the exact update of v.\n");                       exit(EXIT_FAILURE);
  err08: printf("ERROR: LU decomp for computing the determinant failed at the update of R.\n");    exit(EXIT_FAILURE);
}

/* for standard kernels (memory efficient) */
void interpolate1(
  double         *T,
  double         *V,
  const double   *Y,
  const int       N,
  const double   *x,
  const double   *y,
  const double   *w,
  const double   *s,
  const double   *R,
  const double   *t,
  const double   *r,
  const pwsz     sz,
  const pwpm     pm
  ){

  int d,i,j,k,m,n; int D=sz.D,K=sz.K,M=sz.M; int *wi; double *u,*ix,*A,*B,*E,*G,*L,*Q,*W,*wd; char uplo='U'; int *U;
  int info; double bet=pm.bet,lmd=pm.lmd; double val,cc=pm.lmd*SQ(*r/(*s)); int si=sizeof(int),sd=sizeof(double);
  int *ipiv; double *lw; char type=pm.opt&PW_OPT_AFFIN?'a':'s';

  /* allocation */
  W =calloc(D*M,sd); ix=calloc(D*M,sd); ipiv=calloc(D,  sd);
  E =calloc(D*M,sd); u =calloc(D*N,sd); lw  =calloc(D*D,sd);
  /* switch: non-rigid and rigid */
  if(lmd>=1e8||(pm.opt&PW_OPT_NONRG)){for(d=0;d<D;d++)for(n=0;n<N;n++){{u[d+D*n]=Y[d+D*n];}} goto skip;}
  /* non-rigid */
  cc=lmd*SQ(*r/(*s)); invlinear(ix,x,D,M,s,R,t,type,lw,ipiv);
  for(d=0;d<D;d++)for(m=0;m<M;m++) E[d+D*m]=W[m+M*d]=ix[d+D*m]-y[d+D*m];
  if(K){ /* nystrom */
    /* coefficient: W */
    A=calloc(K*K,sd); L =calloc(K,sd); wd=calloc(K*(K+11),sd);
    B=calloc(K*D,sd); U =calloc(N,si);
    Q=calloc(M*K,sd); wi=calloc(M,si);
    gramdecomp(Q,L,wd,wi,y,D,M,K,bet,kernel[pm.G]);
    #pragma omp parallel for private (j) private (m) private (val)
    for(i=0;i<K;i++)for(j=i;j<K;j++){val=(i==j?cc/L[i]:0);for(m=0;m<M;m++){val+=w[m]*Q[m+M*i]*Q[m+M*j];} A[i+K*j]=A[j+K*i]=val;}
    for(d=0;d<D;d++)for(m=0;m<M;m++) E[d+D*m]*=w[m];
    for(k=0;k<K;k++)for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=Q[m+M*k]*E[d+D*m];} B[k+K*d]=val;}
    dposv_(&uplo,&K,&D,A,&K,B,&K,&info); assert(!info);
    for(k=0;k<K;k++)for(m=0;m<M;m++) Q[m+M*k]*=w[m];
    for(m=0;m<M;m++)for(d=0;d<D;d++){val=E[d+D*m];for(k=0;k<K;k++){val-=Q[m+M*k]*B[k+K*d];} W[m+M*d]=val/cc;}
    /* interpolation */
    randperm(U,N);
    for(i=0;i<K;i++)for(j=i;j<K;j++) A[i+K*j]=A[j+K*i]=kernel[pm.G](Y+D*U[i],Y+D*U[j],D,bet)+(i==j?1e-9:0);
    for(k=0;k<K;k++)for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=kernel[pm.G](Y+D*U[k],y+D*m,D,bet)*W[m+M*d];} B[k+K*d]=val;}
    dposv_(&uplo,&K,&D,A,&K,B,&K,&info); assert(!info);
    for(d=0;d<D;d++)for(n=0;n<N;n++){val=0;for(k=0;k<K;k++){val+=kernel[pm.G](Y+D*n,Y+D*U[k],D,bet)*B[k+K*d];} V[d+D*n]=val; u[d+D*n]=val+Y[d+D*n];}
    free(A);free(Q);free(wd);free(U);
    free(B);free(L);free(wi);
  }
  else { /* direct */
    /* coefficient: W */
    G=calloc(M*M,sd);
    #pragma omp parallel for private (j)
    for(i=0;i<M;i++)for(j=i;j<M;j++) G[i+M*j]=G[j+M*i]=kernel[pm.G](y+D*i,y+D*j,D,bet)+(i==j?cc/w[i]:0);
    dposv_(&uplo,&M,&D,G,&M,W,&M,&info); assert(!info);
    /* interpolation */
    #pragma omp parallel for private (d) private (m) private (val)
    for(n=0;n<N;n++)for(d=0;d<D;d++){val=0;for(m=0;m<M;m++){val+=kernel[pm.G](Y+D*n,y+D*m,D,bet)*W[m+M*d];} V[d+D*n]=val; u[d+D*n]=val+Y[d+D*n];}
    free(G);
  }
  skip:
  for(d=0;d<D;d++)for(n=0;n<N;n++){val=0;for(i=0;i<D;i++){val+=R[d+D*i]*u[i+D*n];} T[d+D*n]=(*s)*val+t[d];}
  free(W);free(ix);free(E);free(u);free(ipiv);free(lw);
}

/* for geodesic kernel and ssm (dld++) */
void interpolate2(
  double         *T,
  double         *V,
  const double   *Y,
  const int       N,
  const double   *x,
  const double   *y,
  const double   *w,
  const double   *s,
  const double   *R,
  const double   *t,
  const double   *r,
  const double   *LQ,
  const int      *U,
  const pwsz     sz,
  const pwpm     pm
  ){

  int d,i,k,m,n; int D=sz.D,K=sz.K,M=sz.M; double *u,*ix,*B,*E,*wd;
  const double lmd=pm.lmd; double val; int sd=sizeof(double),si=sizeof(int);
  const double *L=LQ,*Q=LQ+K; double *Qy; int ssm=strlen(pm.fn[COV_LQ]);
  int *ipiv; double *lw; char type=pm.opt&PW_OPT_AFFIN?'a':'s';

  assert(K);

  /* allocation */
  lw=calloc(D*D,sd); ix=calloc(D*M,sd);
  E =calloc(D*M,sd); u =calloc(D*N,sd);
  B =calloc(K*D,sd); Qy=calloc(K*M*(ssm?D:1),sd);
  wd=calloc((ssm?D:1)*M*K+M*D+2*K*K,sd); ipiv=calloc(D,si);

  if(ssm) for(k=0;k<K;k++)for(m=0;m<M;m++)for(d=0;d<D;d++) Qy[d+D*m+D*M*k]=Q[d+D*U[m]+D*N*k];
  else    for(k=0;k<K;k++)for(m=0;m<M;m++) Qy[m+M*k]=Q[U[m]+N*k];
  invlinear(ix,x,D,M,s,R,t,type,lw,ipiv);
  for(m=0;m<M;m++)for(d=0;d<D;d++) E[d+D*m]=ix[d+D*m]-y[d+D*m];

  /* switch: non-rigid and rigid */
  if(lmd>=1e8||(pm.opt&PW_OPT_NONRG)){for(d=0;d<D;d++)for(n=0;n<N;n++){{u[d+D*n]=Y[d+D*n];}} goto skip;}

  lowrank_encode(B,wd,E,w,s,r,L,Qy,D,M,K,lmd,ssm);
  lowrank_decode(u,V,Y,Q,B,D,N,K,ssm);

  skip:
  for(d=0;d<D;d++)for(n=0;n<N;n++){val=0;for(i=0;i<D;i++){val+=R[d+D*i]*u[i+D*n];} T[d+D*n]=(*s)*val+t[d];}

  free(B); free(Qy); free(wd); free(ipiv);
  free(u); free(ix); free(E);  free(lw);
}

void interpolate_x(
    double       *x,
    const double *y,
    const double *X,
    int           D,
    int           M,
    int           N,
    const double  r,
    pwpm          pm
  ){
  int d,k,m,K=30; double e,val; double *p,*P; int *T,*q,*Q;
  int me=-1; int si=sizeof(int),sd=sizeof(double);

  K=K<M?K:M;

  P=calloc(M*(K+1),sd); T=kdtree_build(X,D,N);
  Q=calloc(M*(K+1),si);

  #pragma omp parallel for private (q) private (p) private (d) private (e) private (k)
  for(m=0;m<M;m++){q=Q+(K+1)*m;p=P+(K+1)*m; knnsearch(q,K,2*r,y+D*m,me,X,T,D,N);
    if(*q){*p=0;for(k=1;k<=*q;k++){p[k]=gauss(y+D*m,X+D*q[k],D,r);*p+=p[k];}}
    else{*p=1.0;p[1]=1.0;*q=1;nnsearch(q+1,&e,y+D*m,X,T,D,N);}
    for(d=0;d<D;d++){val=0;for(k=1;k<=*q;k++){val+=p[k]*X[d+D*q[k]];} x[d+D*m]=val/(*p);}
  }
  free(P);free(Q);free(T);

  return;
}

