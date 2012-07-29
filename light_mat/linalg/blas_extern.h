/**
 * @file blas_extern.h
 *
 * External function declarations for BLAS
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BLAS_EXTERN_H_
#define LIGHTMAT_BLAS_EXTERN_H_

#include <light_mat/config/config.h>


/************************************************
 *
 *  function name macros
 *
 ************************************************/

// BLAS Level 1

#define LMAT_SASUM   sasum
#define LMAT_SAXPY   saxpy
#define LMAT_SDOT    sdot
#define LMAT_SNRM2   snrm2
#define LMAT_SROT    srot

#define LMAT_DASUM   dasum
#define LMAT_DAXPY   daxpy
#define LMAT_DDOT    ddot
#define LMAT_DNRM2   dnrm2
#define LMAT_DROT    drot

// BLAS Level 2

#define LMAT_SGEMV	sgemv
#define LMAT_SGER	sger
#define LMAT_SSYMV	ssymv
#define LMAT_STRMV	strmv
#define LMAT_STRSV	strsv

#define LMAT_DGEMV	dgemv
#define LMAT_DGER	dger
#define LMAT_DSYMV	dsymv
#define LMAT_DTRMV	dtrmv
#define LMAT_DTRSV 	dtrsv


// BLAS Level 3

#define LMAT_SGEMM	sgemm
#define LMAT_SSYMM	ssymm
#define LMAT_STRMM 	strmm
#define LMAT_STRSM 	strsm

#define LMAT_DGEMM	dgemm
#define LMAT_DSYMM	dsymm
#define LMAT_DTRMM	dtrmm
#define LMAT_DTRSM	dtrsm

typedef int lmat_blas_int;

extern "C"
{
	/************************************************
	 *
	 *  BLAS Level 1
	 *
	 ************************************************/

	float   LMAT_SASUM (const lmat_blas_int *n, const float *x, const lmat_blas_int *incx);

	void    LMAT_SAXPY (const lmat_blas_int *n,
						const float *alpha, const float *x, const lmat_blas_int *incx,
						float *y, const lmat_blas_int *incy);

	float   LMAT_SDOT  (const lmat_blas_int *n,
						const float *x, const lmat_blas_int *incx,
						const float *y, const lmat_blas_int *incy);

	float   LMAT_SNRM2 (const lmat_blas_int *n,
						const float *x, const lmat_blas_int *incx);

	void    LMAT_SROT  (const lmat_blas_int *n,
						float *x, const lmat_blas_int *incx,
						float *y, const lmat_blas_int *incy, const float *c, const float *s);

	double  LMAT_DASUM (const lmat_blas_int *n,
						const double *x, const lmat_blas_int *incx);

	void    LMAT_DAXPY (const lmat_blas_int *n,
						const double *alpha, const double *x, const lmat_blas_int *incx,
						double *y, const lmat_blas_int *incy);

	double  LMAT_DDOT  (const lmat_blas_int *n,
						const double *x, const lmat_blas_int *incx,
						const double *y, const lmat_blas_int *incy);

	double  LMAT_DNRM2 (const lmat_blas_int *n, const double *x, const lmat_blas_int *incx);

	void    LMAT_DROT  (const lmat_blas_int *n,
						double *x, const lmat_blas_int *incx,
						double *y, const lmat_blas_int *incy, const double *c, const double *s);


	/************************************************
	 *
	 *  BLAS Level 2
	 *
	 ************************************************/

	void LMAT_SGEMV (const char *trans, const lmat_blas_int *m, const lmat_blas_int *n,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
					 const float *x, const lmat_blas_int *incx,
	           	     const float *beta, float *y, const lmat_blas_int *incy);

	void LMAT_SGER  (const lmat_blas_int *m, const lmat_blas_int *n,
					 const float *alpha, const float *x, const lmat_blas_int *incx,
					 const float *y, const lmat_blas_int *incy,
					 float *a, const lmat_blas_int *lda);

	void LMAT_SSYMV (const char *uplo, const lmat_blas_int *n,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
	           	     const float *x, const lmat_blas_int *incx, const float *beta,
	           	     float *y, const lmat_blas_int *incy);

	void LMAT_STRMV (const char *uplo, const char *transa,
					 const char *diag, const lmat_blas_int *n,
	           	     const float *a, const lmat_blas_int *lda,
	           	     float *b, const lmat_blas_int *incx);

	void LMAT_STRSV (const char *uplo, const char *trans,
					 const char *diag, const lmat_blas_int *n,
	           	     const float *a, const lmat_blas_int *lda,
	           	     float *x, const lmat_blas_int *incx);

	void LMAT_DGEMV (const char *trans, const lmat_blas_int *m, const lmat_blas_int *n,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
					 const double *x, const lmat_blas_int *incx,
	           	     const double *beta, double *y, const lmat_blas_int *incy);

	void LMAT_DGER  (const lmat_blas_int *m, const lmat_blas_int *n,
					 const double *alpha, const double *x, const lmat_blas_int *incx,
	          	     const double *y, const lmat_blas_int *incy,
	          	     double *a, const lmat_blas_int *lda);

	void LMAT_DSYMV (const char *uplo, const lmat_blas_int *n,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
	           	     const double *x, const lmat_blas_int *incx,
	           	     const double *beta, double *y, const lmat_blas_int *incy);

	void LMAT_DTRMV (const char *uplo, const char *transa,
					 const char *diag, const lmat_blas_int *n,
	           	     const double *a, const lmat_blas_int *lda,
	           	     double *b, const lmat_blas_int *incx);

	void LMAT_DTRSV (const char *uplo, const char *trans,
					 const char *diag, const lmat_blas_int *n,
					 const double *a, const lmat_blas_int *lda,
					 double *x, const lmat_blas_int *incx);


	/************************************************
	 *
	 *  BLAS Level 3
	 *
	 ************************************************/

	void LMAT_SGEMM (const char *transa, const char *transb,
					 const lmat_blas_int *m, const lmat_blas_int *n, const lmat_blas_int *k,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
					 const float *b, const lmat_blas_int *ldb,
					 const float *beta, float *c, const lmat_blas_int *ldc);

	void LMAT_SSYMM (const char *side, const char *uplo,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
					 const float *b, const lmat_blas_int *ldb,
					 const float *beta, float *c, const lmat_blas_int *ldc);

	void LMAT_STRMM (const char *side, const char *uplo,
					 const char *transa, const char *diag,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
					 float *b, const lmat_blas_int *ldb);

	void LMAT_STRSM (const char *side, const char *uplo,
					 const char *transa, const char *diag,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const float *alpha, const float *a, const lmat_blas_int *lda,
					 float *b, const lmat_blas_int *ldb);

	void LMAT_DGEMM (const char *transa, const char *transb,
					 const lmat_blas_int *m, const lmat_blas_int *n, const lmat_blas_int *k,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
					 const double *b, const lmat_blas_int *ldb,
					 const double *beta, double *c, const lmat_blas_int *ldc);

	void LMAT_DSYMM (const char *side, const char *uplo,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
					 const double *b, const lmat_blas_int *ldb,
					 const double *beta, double *c, const lmat_blas_int *ldc);

	void LMAT_DTRMM (const char *side, const char *uplo,
					 const char *transa, const char *diag,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
					 double *b, const lmat_blas_int *ldb);

	void LMAT_DTRSM (const char *side, const char *uplo,
					 const char *transa, const char *diag,
					 const lmat_blas_int *m, const lmat_blas_int *n,
					 const double *alpha, const double *a, const lmat_blas_int *lda,
					 double *b, const lmat_blas_int *ldb);
}

#endif /* BLAS_EXTERN_H_ */
