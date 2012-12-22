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

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/linalg/blas_extern.h>
#include "internal/linalg_aux.h"

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
