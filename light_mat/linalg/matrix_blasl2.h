/**
 * @file matrix_blasl2.h
 *
 * BLAS Level 2 on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_BLASL2_H_
#define LIGHTMAT_MATRIX_BLASL2_H_

#include "bits/matrix_blasl2_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  gemv_n
	 *
	 ********************************************/

	template<typename T, class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void check_gemv_n_args(
			const IMatrixXpr<MatA, T>& A,
			const IMatrixXpr<VecX, T>& x,
			IDenseMatrix<VecY, T>& y)
	{
		check_arg(is_column(x) && x.nrows() == A.ncolumns(),
				"blas::gemv_n: x should be a column vector of length n");

		check_arg(is_column(y) && y.nrows() == A.nrows(),
				"blas::gemv_n: y should be a column vector of length m");
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			IDenseMatrix<VecY, float>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(A, x, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			IDenseMatrix<VecY, double>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(A, x, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(const float alpha,
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			IDenseMatrix<VecY, float>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(alpha, A, x, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(const double alpha,
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			IDenseMatrix<VecY, double>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(alpha, A, x, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(A, x, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(A, x, beta, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(const float alpha,
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(alpha, A, x, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_n(const double alpha,
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gemv_n_args(A, x, y);

		typedef typename lmat::detail::gemv_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_n(alpha, A, x, beta, y);
	}


	/********************************************
	 *
	 *  gemv_t
	 *
	 ********************************************/

	template<typename T, class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void check_gemv_t_args(
			const IMatrixXpr<MatA, T>& A,
			const IMatrixXpr<VecX, T>& x,
			IDenseMatrix<VecY, T>& y)
	{
		check_arg(is_column(x) && x.nrows() == A.nrows(),
				"blas::gemv_t: x should be a column vector of length m");

		check_arg(is_column(y) && y.nrows() == A.ncolumns(),
				"blas::gemv_t: y should be a column vector of length n");
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			IDenseMatrix<VecY, float>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(A, x, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			IDenseMatrix<VecY, double>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(A, x, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(const float alpha,
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			IDenseMatrix<VecY, float>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(alpha, A, x, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(const double alpha,
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			IDenseMatrix<VecY, double>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(alpha, A, x, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(A, x, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(A, x, beta, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(const float alpha,
			const IMatrixXpr<MatA, float>& A,
			const IMatrixXpr<VecX, float>& x,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(alpha, A, x, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gemv_t(const double alpha,
			const IMatrixXpr<MatA, double>& A,
			const IMatrixXpr<VecX, double>& x,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gemv_t_args(A, x, y);

		typedef typename lmat::detail::gemv_t_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval_t(alpha, A, x, beta, y);
	}


}

#endif /* MATRIX_BLASL2_H_ */
