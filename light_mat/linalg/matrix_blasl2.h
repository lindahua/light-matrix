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
		intern_t::eval(A, x, y);
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
		intern_t::eval(A, x, y);
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
		intern_t::eval(alpha, A, x, y);
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
		intern_t::eval(alpha, A, x, y);
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
		intern_t::eval(A, x, beta, y);
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
		intern_t::eval(A, x, beta, y);
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
		intern_t::eval(alpha, A, x, beta, y);
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
		intern_t::eval(alpha, A, x, beta, y);
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
		intern_t::eval(A, x, y);
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
		intern_t::eval(A, x, y);
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
		intern_t::eval(alpha, A, x, y);
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
		intern_t::eval(alpha, A, x, y);
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
		intern_t::eval(A, x, beta, y);
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
		intern_t::eval(A, x, beta, y);
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
		intern_t::eval(alpha, A, x, beta, y);
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
		intern_t::eval(alpha, A, x, beta, y);
	}


	/********************************************
	 *
	 *  gevm_n
	 *
	 ********************************************/

	template<typename T, class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void check_gevm_n_args(
			const IMatrixXpr<VecX, T>& x,
			const IMatrixXpr<MatA, T>& A,
			IDenseMatrix<VecY, T>& y)
	{
		check_arg(is_row(x) && x.ncolumns() == A.nrows(),
				"blas::gemv_n: x should be a row vector of length m");

		check_arg(is_row(y) && y.ncolumns() == A.ncolumns(),
				"blas::gemv_n: y should be a column vector of length n");
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(
			const IMatrixXpr<VecX, float>& x,
			const IMatrixXpr<MatA, float>& A,
			IDenseMatrix<VecY, float>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(x, A, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(
			const IMatrixXpr<VecX, double>& x,
			const IMatrixXpr<MatA, double>& A,
			IDenseMatrix<VecY, double>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(x, A, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(const float alpha,
			const IMatrixXpr<VecX, float>& x,
			const IMatrixXpr<MatA, float>& A,
			IDenseMatrix<VecY, float>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(alpha, x, A, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(const double alpha,
			const IMatrixXpr<VecX, double>& x,
			const IMatrixXpr<MatA, double>& A,
			IDenseMatrix<VecY, double>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(alpha, x, A, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(
			const IMatrixXpr<VecX, float>& x,
			const IMatrixXpr<MatA, float>& A,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(x, A, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(
			const IMatrixXpr<VecX, double>& x,
			const IMatrixXpr<MatA, double>& A,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(x, A, beta, y);
	}


	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(const float alpha,
			const IMatrixXpr<VecX, float>& x,
			const IMatrixXpr<MatA, float>& A,
			const float beta, IDenseMatrix<VecY, float>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(alpha, x, A, beta, y);
	}

	template<class MatA, class VecX, class VecY>
	LMAT_ENSURE_INLINE
	void gevm_n(const double alpha,
			const IMatrixXpr<VecX, double>& x,
			const IMatrixXpr<MatA, double>& A,
			const double beta, IDenseMatrix<VecY, double>& y)
	{
		check_gevm_n_args(x, A, y);

		typedef typename lmat::detail::gevm_n_internal_map<MatA, VecX, VecY>::type intern_t;
		intern_t::eval(alpha, x, A, beta, y);
	}


}

#endif /* MATRIX_BLASL2_H_ */
