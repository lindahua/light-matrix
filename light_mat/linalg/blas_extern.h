/**
 * @file blas_extern.h
 *
 * External import of BLAS functions
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BLAS_EXTERN_H_
#define LIGHTMAT_BLAS_EXTERN_H_

#include <light_mat/config/config.h>

#ifdef LMAT_BLAS_ILP64
typedef long long blas_int;
#else
typedef int blas_int;
#endif

#define LMAT_BLAS_NAME(name) name

extern "C"
{
	// BLAS LEVEL 1

	float   LMAT_BLAS_NAME(sasum)(const blas_int *n, const float *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(saxpy)(const blas_int *n, const float *alpha, const float *x, const blas_int *incx, float *y, const blas_int *incy);
	void    LMAT_BLAS_NAME(scopy)(const blas_int *n, const float *x, const blas_int *incx, float *y, const blas_int *incy);
	float   LMAT_BLAS_NAME(sdot) (const blas_int *n, const float *x, const blas_int *incx, const float *y, const blas_int *incy);
	float   LMAT_BLAS_NAME(snrm2)(const blas_int *n, const float *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(srot) (const blas_int *n, float *x, const blas_int *incx, float *y, const blas_int *incy, const float *c, const float *s);
	void    LMAT_BLAS_NAME(sscal)(const blas_int *n, const float *a, float *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(sswap)(const blas_int *n, float *x, const blas_int *incx, float *y, const blas_int *incy);

	double  LMAT_BLAS_NAME(dasum)(const blas_int *n, const double *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(daxpy)(const blas_int *n, const double *alpha, const double *x, const blas_int *incx, double *y, const blas_int *incy);
	void    LMAT_BLAS_NAME(dcopy)(const blas_int *n, const double *x, const blas_int *incx, double *y, const blas_int *incy);
	double  LMAT_BLAS_NAME(ddot) (const blas_int *n, const double *x, const blas_int *incx, const double *y, const blas_int *incy);
	double  LMAT_BLAS_NAME(dnrm2)(const blas_int *n, const double *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(drot) (const blas_int *n, double *x, const blas_int *incx, double *y, const blas_int *incy, const double *c, const double *s);
	void    LMAT_BLAS_NAME(dscal)(const blas_int *n, const double *a, double *x, const blas_int *incx);
	void    LMAT_BLAS_NAME(dswap)(const blas_int *n, double *x, const blas_int *incx, double *y, const blas_int *incy);

	// BLAS LEVEL 2

	void LMAT_BLAS_NAME(sgemv)(const char *trans, const blas_int *m, const blas_int *n, const float *alpha,
	           const float *a, const blas_int *lda, const float *x, const blas_int *incx,
	           const float *beta, float *y, const blas_int *incy);
	void LMAT_BLAS_NAME(sger)(const blas_int *m, const blas_int *n, const float *alpha, const float *x, const blas_int *incx,
	          const float *y, const blas_int *incy, float *a, const blas_int *lda);
	void LMAT_BLAS_NAME(ssymv)(const char *uplo, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           const float *x, const blas_int *incx, const float *beta, float *y, const blas_int *incy);
	void LMAT_BLAS_NAME(ssyr)(const char *uplo, const blas_int *n, const float *alpha, const float *x, const blas_int *incx,
	          float *a, const blas_int *lda);
	void LMAT_BLAS_NAME(strmv)(const char *uplo, const char *transa, const char *diag, const blas_int *n, const float *a,
	           const blas_int *lda, float *b, const blas_int *incx);
	void LMAT_BLAS_NAME(strsv)(const char *uplo, const char *trans, const char *diag, const blas_int *n,
	           const float *a, const blas_int *lda, float *x, const blas_int *incx);

	void LMAT_BLAS_NAME(dgemv)(const char *trans, const blas_int *m, const blas_int *n, const double *alpha,
	           const double *a, const blas_int *lda, const double *x, const blas_int *incx,
	           const double *beta, double *y, const blas_int *incy);
	void LMAT_BLAS_NAME(dger)(const blas_int *m, const blas_int *n, const double *alpha, const double *x, const blas_int *incx,
	          const double *y, const blas_int *incy, double *a, const blas_int *lda);
	void LMAT_BLAS_NAME(dsymv)(const char *uplo, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           const double *x, const blas_int *incx, const double *beta, double *y, const blas_int *incy);
	void LMAT_BLAS_NAME(dsyr)(const char *uplo, const blas_int *n, const double *alpha, const double *x, const blas_int *incx,
	          double *a, const blas_int *lda);
	void LMAT_BLAS_NAME(dtrmv)(const char *uplo, const char *transa, const char *diag, const blas_int *n,
	           const double *a, const blas_int *lda, double *b, const blas_int *incx);
	void LMAT_BLAS_NAME(dtrsv)(const char *uplo, const char *trans, const char *diag, const blas_int *n,
	           const double *a, const blas_int *lda, double *x, const blas_int *incx);

	// BLAS LEVEL 3

	void LMAT_BLAS_NAME(sgemm)(const char *transa, const char *transb, const blas_int *m, const blas_int *n, const blas_int *k,
	           const float *alpha, const float *a, const blas_int *lda, const float *b, const blas_int *ldb,
	           const float *beta, float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(ssymm)(const char *side, const char *uplo, const blas_int *m, const blas_int *n,
	           const float *alpha, const float *a, const blas_int *lda, const float *b, const blas_int *ldb,
	           const float *beta, float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(ssyrk)(const char *uplo, const char *trans, const blas_int *n, const blas_int *k,
	           const float *alpha, const float *a, const blas_int *lda, const float *beta,
	           float *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(strmm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           float *b, const blas_int *ldb);
	void LMAT_BLAS_NAME(strsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const float *alpha, const float *a, const blas_int *lda,
	           float *b, const blas_int *ldb);

	void LMAT_BLAS_NAME(dgemm)(const char *transa, const char *transb, const blas_int *m, const blas_int *n, const blas_int *k,
	           const double *alpha, const double *a, const blas_int *lda, const double *b, const blas_int *ldb,
	           const double *beta, double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dsymm)(const char *side, const char *uplo, const blas_int *m, const blas_int *n,
	           const double *alpha, const double *a, const blas_int *lda, const double *b, const blas_int *ldb,
	           const double *beta, double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dsyrk)(const char *uplo, const char *trans, const blas_int *n, const blas_int *k,
	           const double *alpha, const double *a, const blas_int *lda, const double *beta,
	           double *c, const blas_int *ldc);
	void LMAT_BLAS_NAME(dtrmm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           double *b, const blas_int *ldb);
	void LMAT_BLAS_NAME(dtrsm)(const char *side, const char *uplo, const char *transa, const char *diag,
	           const blas_int *m, const blas_int *n, const double *alpha, const double *a, const blas_int *lda,
	           double *b, const blas_int *ldb);

}




#endif 
