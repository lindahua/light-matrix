/**
 * @file blas_l2.h
 *
 * @brief BLAS Level 2 functions
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_BLAS_L2_H_
#define LIGHTMAT_BLAS_L2_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/linalg/blas_extern.h>
#include "internal/linalg_aux.h"

namespace lmat { namespace blas {

	// gemv

	namespace internal
	{
		template<class A, class X, class Y>
		LMAT_ENSURE_INLINE
		inline bool gemv_check_dims(const A& a, const X& x, const Y& y, char trans)
		{
			if (trans == 'N' || trans == 'n')
			{
				return a.ncolumns() == x.nelems() && a.nrows() == y.nelems();
			}
			else
			{
				return a.nrows() == x.nelems() && a.ncolumns() == y.nelems();
			}
		}
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			float alpha, const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x,
			float beta, IRegularMatrix<Y, float>& y, char trans='N')
	{
		LMAT_CHECK_DIMS( internal::gemv_check_dims(a, x, y, trans) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int m = a.nrows();
		blas_int n = a.ncolumns();

		LMAT_BLAS_NAME(sgemv)(&trans, &m, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			double alpha, const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x,
			double beta, IRegularMatrix<Y, double>& y, char trans='N')
	{
		LMAT_CHECK_DIMS( internal::gemv_check_dims(a, x, y, trans) )

		blas_int lda = (blas_int)a.col_stride();
		blas_int incx = (blas_int)lmat::internal::get_vector_intv(x);
		blas_int incy = (blas_int)lmat::internal::get_vector_intv(y);

		blas_int m = a.nrows();
		blas_int n = a.ncolumns();

		LMAT_BLAS_NAME(dgemv)(&trans, &m, &n, &alpha, a.ptr_data(), &lda, x.ptr_data(), &incx, &beta, y.ptr_data(), &incy);
	}


	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			const IRegularMatrix<A, float>& a, const IRegularMatrix<X, float>& x, IRegularMatrix<Y, float>& y, char trans='N')
	{
		gemv(1.0f, a, x, 0.0f, y, trans);
	}

	template<class A, class X, class Y>
	LMAT_ENSURE_INLINE
	inline void gemv(
			const IRegularMatrix<A, double>& a, const IRegularMatrix<X, double>& x, IRegularMatrix<Y, double>& y, char trans='N')
	{
		gemv(1.0, a, x, 0.0, y, trans);
	}


} }


#endif /* BLAS_L2_H_ */
