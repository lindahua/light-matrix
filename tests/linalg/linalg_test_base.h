/**
 * @file linalg_test_base.h
 *
 * The testing facilities for linear algebra unit testing
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LINALG_TEST_BASE_H_
#define LIGHTMAT_LINALG_TEST_BASE_H_

#include "../test_base.h"
#include "../multimat_supp.h"
#include <light_mat/math/math_base.h>

namespace lmat { namespace test {

	template<typename T> struct blas_default_tol;

	template<> struct blas_default_tol<float>
	{
		static float get()
		{
			return 1.0e-5f;
		}
	};

	template<> struct blas_default_tol<double>
	{
		static double get()
		{
			return 1.0e-12;
		}
	};

	template<class Mat, typename T>
	index_t safe_get_vec_intv(const IRegularMatrix<Mat, T>& mat)
	{
		index_t m = mat.nrows();
		index_t n = mat.ncolumns();
		index_t rs = mat.row_stride();
		index_t cs = mat.col_stride();

		if (cs == m && rs == 1)
		{
			return 1;
		}
		else if (n == 1)
		{
			return rs;
		}
		else if (m == 1)
		{
			return cs;
		}
		else
		{
			throw invalid_argument("safe_get_vec_intv:not a proper vector.");
		}
	}


	/********************************************
	 *
	 *  emulation of BLAS functions
	 *
	 ********************************************/

	// Level 1

	template<typename T>
	T safe_asum(index_t len, const T *x, index_t incx)
	{
		double s(0);
		for (index_t i = 0; i < len; ++i)
		{
			s += math::abs(x[i * incx]);
		}
		return static_cast<T>(s);
	}

	template<typename T>
	T safe_dot(index_t len, const T *x, index_t incx, const T *y, index_t incy)
	{
		double s(0);
		for (index_t i = 0; i < len; ++i)
		{
			s += x[i * incx] * y[i * incy];
		}
		return static_cast<T>(s);
	}

	template<typename T>
	T safe_nrm2(index_t len, const T *x, index_t incx)
	{
		return math::sqrt(safe_dot(len, x, incx, x, incx));
	}

	template<typename T>
	void safe_scal(index_t len, T a, const T *x, index_t incx, T *r)
	{
		for (index_t i = 0; i < len; ++i)
		{
			r[i] = a * x[i * incx];
		}
	}

	template<typename T>
	void safe_axpy(index_t len, T a, const T *x, index_t incx, const T *y, index_t incy, T *r)
	{
		for (index_t i = 0; i < len; ++i)
		{
			r[i] = a * x[i * incx] + y[i * incy];
		}
	}

	template<typename T>
	void safe_rot(index_t len, T c, T s, const T *x, index_t incx, const T *y, index_t incy, T *rx, T *ry)
	{
		for (index_t i = 0; i < len; ++i)
		{
			T xi = x[i * incx];
			T yi = y[i * incy];

			rx[i] = c * xi + s * yi;
			ry[i] = c * yi - s * xi;
		}
	}


	// Level 2 & 3

	template<typename T, class A, class X, class Y, class R>
	void safe_mv(const T& alpha, const A& a, char trans, const X& x, const T& beta, const Y& y, R& r)
	{
		index_t m = a.nrows();
		index_t n = a.ncolumns();

		if (trans == 'n' || trans == 'N')
		{
			// (m x n) . (n) --> (m)

			for (index_t i = 0; i < m ; ++i)
			{
				T v = 0;
				for (index_t j = 0; j < n; ++j)
					v += a(i, j) * x[j];

				r[i] = alpha * v + beta * y[i];
			}
		}
		else
		{
			// (n x m) . (m) --> (n)

			for (index_t j = 0; j < n; ++j)
			{
				T v = 0;
				for (index_t i = 0; i < m; ++i)
					v += a(i, j) * x[i];

				r[j] = alpha * v + beta * y[j];
			}
		}
	}

} }

#endif

