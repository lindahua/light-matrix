/**
 * @file blas_l1.h
 *
 * BLAS level-1 functions on matrices
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BLAS_L1_H_
#define LIGHTMAT_BLAS_L1_H_

#include "internal/linalg_aux.h"

extern "C"
{
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
}


namespace lmat { namespace blas {

	// asum

	template<class X>
	LMAT_ENSURE_INLINE
	inline float asum(const IRegularMatrix<X, float>& x)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		return LMAT_BLAS_NAME(sasum)(&n, x.ptr_data(), &incx);
	}

	template<class X>
	LMAT_ENSURE_INLINE
	inline double asum(const IRegularMatrix<X, double>& x)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		return LMAT_BLAS_NAME(dasum)(&n, x.ptr_data(), &incx);
	}


	// axpy

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline void axpy(float a, const IRegularMatrix<X, float>& x, IRegularMatrix<Y, float>& y)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		LMAT_BLAS_NAME(saxpy)(&n, &a, x.ptr_data(), &incx, y.ptr_data(), &incy);
	}

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline void axpy(double a, const IRegularMatrix<X, double>& x, IRegularMatrix<Y, double>& y)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		LMAT_BLAS_NAME(daxpy)(&n, &a, x.ptr_data(), &incx, y.ptr_data(), &incy);
	}


	// nrm2

	template<class X>
	LMAT_ENSURE_INLINE
	inline float nrm2(const IRegularMatrix<X, float>& x)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		return LMAT_BLAS_NAME(snrm2)(&n, x.ptr_data(), &incx);
	}

	template<class X>
	LMAT_ENSURE_INLINE
	inline double nrm2(const IRegularMatrix<X, double>& x)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		return LMAT_BLAS_NAME(dnrm2)(&n, x.ptr_data(), &incx);
	}

	// dot

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline float dot(const IRegularMatrix<X, float>& x, const IRegularMatrix<Y, float>& y)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		return LMAT_BLAS_NAME(sdot)(&n, x.ptr_data(), &incx, y.ptr_data(), &incy);
	}

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline double dot(const IRegularMatrix<X, double>& x, const IRegularMatrix<Y, double>& y)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		return LMAT_BLAS_NAME(ddot)(&n, x.ptr_data(), &incx, y.ptr_data(), &incy);
	}

	// rot

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline void rot(IRegularMatrix<X, float>& x, IRegularMatrix<Y, float>& y, float c, float s)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		LMAT_BLAS_NAME(srot)(&n, x.ptr_data(), &incx, y.ptr_data(), &incy, &c, &s);
	}

	template<class X, class Y>
	LMAT_ENSURE_INLINE
	inline void rot(IRegularMatrix<X, double>& x, IRegularMatrix<Y, double>& y, double c, double s)
	{
		blas_int n = (blas_int)(x.nelems());
		LMAT_CHECK_DIMS( x.nelems() == y.nelems() );

		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));
		blas_int incy = (blas_int)(lmat::internal::get_vector_intv(y));

		LMAT_BLAS_NAME(drot)(&n, x.ptr_data(), &incx, y.ptr_data(), &incy, &c, &s);
	}

	// scale

	template<class X>
	LMAT_ENSURE_INLINE
	inline void scal(IRegularMatrix<X, float>& x, float a)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		LMAT_BLAS_NAME(sscal)(&n, &a, x.ptr_data(), &incx);
	}

	template<class X>
	LMAT_ENSURE_INLINE
	inline void scal(IRegularMatrix<X, double>& x, double a)
	{
		blas_int n = (blas_int)(x.nelems());
		blas_int incx = (blas_int)(lmat::internal::get_vector_intv(x));

		LMAT_BLAS_NAME(dscal)(&n, &a, x.ptr_data(), &incx);
	}


} }

#endif 
