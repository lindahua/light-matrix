/**
 * @file matrix_blasl2_internal.h
 *
 * Internal implementation of BLAS Level 2
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL2_INTERNAL_H_
#define LIGHTMAT_MATRIX_BLASL2_INTERNAL_H_

#include <light_mat/linalg/linalg_base.h>
#include <light_mat/linalg/blas_extern.h>

#include "matrix_blas_proxy.h"

namespace lmat { namespace detail {

	/********************************************
	 *
	 *  GEMV
	 *
	 ********************************************/

	template<typename T, int M, int N> struct gemv_internal;

	template<int M, int N>
	struct gemv_internal<double, M, N>
	{
		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_n(
				const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				IDenseMatrix<VecY, double>& y)
		{
			eval_n(1.0f, A, x, 0.0f, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_n(
				const double alpha, const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				IDenseMatrix<VecY, double>& y)
		{
			eval_n(alpha, A, x, 0.0f, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_n(
				const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				const double beta, IDenseMatrix<VecY, double>& y)
		{
			eval_n(1.0f, A, x, beta, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_n(
				const double alpha, const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				const double beta, IDenseMatrix<VecY, double>& y)
		{
			typename lmat::detail::blas_mat_proxy<MatA>::type A_(A.derived());
			typename lmat::detail::blas_col_proxy<VecX>::type x_(x.derived());

			lmat_blas_int m = A_.nrows();
			lmat_blas_int n = A_.ncolumns();
			lmat_blas_int lda = A_.lead_dim();
			lmat_blas_int incx = x_.inc();
			lmat_blas_int incy = 1;

			LMAT_DGEMV("N", &m, &n, &alpha, A_.pdata(), lda, x_.pdata(), &incx, &beta, y.ptr_data(), &incy);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_t(
				const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				IDenseMatrix<VecY, double>& y)
		{
			eval_t(1.0f, A, x, 0.0f, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_t(
				const double alpha, const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				IDenseMatrix<VecY, double>& y)
		{
			eval_t(alpha, A, x, 0.0f, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_t(
				const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				const double beta, IDenseMatrix<VecY, double>& y)
		{
			eval_t(1.0f, A, x, beta, y);
		}

		template<class MatA, class VecX, class VecY>
		LMAT_ENSURE_INLINE
		static void eval_t(
				const double alpha, const IMatrixXpr<MatA, double>& A,
				const IMatrixXpr<VecX, double>& x,
				const double beta, IDenseMatrix<VecY, double>& y)
		{
			typename lmat::detail::blas_mat_proxy<MatA>::type A_(A.derived());
			typename lmat::detail::blas_col_proxy<VecX>::type x_(x.derived());

			lmat_blas_int m = A_.nrows();
			lmat_blas_int n = A_.ncolumns();
			lmat_blas_int lda = A_.lead_dim();
			lmat_blas_int incx = x_.inc();
			lmat_blas_int incy = 1;

			LMAT_DGEMV("T", &m, &n, &alpha, A_.pdata(), lda, x_.pdata(), &incx, &beta, y.ptr_data(), &incy);
		}
	};


	template<class MatA, class VecX, class VecY>
	struct gemv_n_internal_map
	{
		typedef typename matrix_traits<MatA>::value_type T;
		static const int M = binary_ctdim<ct_rows<MatA>::value, ct_rows<VecY>::value>::value;
		static const int N = binary_ctdim<ct_cols<MatA>::value, ct_rows<VecX>::value>::value;

		typedef gemv_internal<T, M, N> type;
	};


} }

#endif /* MATRIX_BLASL2_INTERNAL_H_ */
