int dsyev_ (char *jobz,char *uplo,int *n,double *a,int *lda,double *w,double *work,int *lwork,int *info);
int dposv_ (char *uplo,int *n,int *nrhs,double *A,int *lda,double *B,int *ldb,int *info);
int dpotrs_(char *uplo,int *n,int *nrhs,double *A,int *lda,double *B,int *ldb,int *info);
int dpotrf_(char *uplo,int *n,double *A,int *lda,int *info);
int dpotri_(char *uplo,int *n,double *A,int *lda,int *info);
int dgetrf_(int *m,int *n,double *A,int *lda,int *ipiv,int *info);

