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
	 *  prepare specific type of matrices
	 *
	 ********************************************/

	template<typename T>
	T randunif(const T& a, const T& b)
	{
		double u = double(std::rand()) / RAND_MAX;
		return static_cast<T>(a + (b - a) * u);
	}


	template<typename T, class Mat>
	void fill_rand_sym(IRegularMatrix<Mat, T>& mat)
	{
		index_t n = mat.nrows();

		for (index_t i = 1; i < n; ++i)
		{
			for (index_t j = 0; j < i; ++j)
			{
				T v = randunif<T>(T(-1.0), T(1.0));
				mat(i, j) = v;
				mat(j, i) = v;
			}
		}

		for (index_t i = 0; i < n; ++i)
		{
			mat(i, i) = randunif<T>(T(-1.0), T(1.0));
		}
	}


	template<typename T, class Mat>
	void fill_rand_tri(IRegularMatrix<Mat, T>& mat, char uplo, bool uni=false)
	{
		index_t n = mat.nrows();

		if (uplo == 'L' || uplo == 'l')
		{
			for (index_t i = 1; i < n; ++i)
			{
				for (index_t j = 0; j < i; ++j)
				{
					T v = randunif<T>(T(-0.5), T(0.5));
					mat(i, j) = v;
				}
			}
		}
		else
		{
			for (index_t i = 0; i < n-1; ++i)
			{
				for (index_t j = i+1; j < n; ++j)
				{
					T v = randunif<T>(T(-0.5), T(0.5));
					mat(i, j) = v;
				}
			}
		}

		if (uni)
		{
			for (index_t i = 0; i < n; ++i)
				mat(i, i) = T(1.0);
		}
		else
		{
			for (index_t i = 0; i < n; ++i)
				mat(i, i) = randunif<T>(T(1.0), T(3.0));
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


	// Level 2

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
				double v = 0;
				for (index_t j = 0; j < n; ++j)
					v += a(i, j) * x[j];

				r[i] = (T)(alpha * v + beta * y[i]);
			}
		}
		else
		{
			// (n x m) . (m) --> (n)

			for (index_t j = 0; j < n; ++j)
			{
				double v = 0;
				for (index_t i = 0; i < m; ++i)
					v += a(i, j) * x[i];

				r[j] = (T)(alpha * v + beta * y[j]);
			}
		}
	}


	// Level 3

	template<typename T, class A, class B, class C, class R>
	void safe_mm(const T& alpha, const A& a, char transa, const B& b, char transb,
			const T& beta, const C& c, R& r)
	{
		bool ta = !(transa == 'n' || transa == 'N');
		bool tb = !(transb == 'n' || transb == 'N');

		if (!ta && !tb)
		{
			index_t m = a.nrows();
			index_t k = a.ncolumns();
			index_t n = b.ncolumns();

			for (index_t j = 0; j < n; ++j)
			{
				for (index_t i = 0; i < m; ++i)
				{
					double v = 0;
					for (index_t l = 0; l < k; ++l) v += a(i, l) * b(l, j);
					r(i, j) = (T)(alpha * v + beta * c(i, j));
				}
			}
		}
		else if (!ta && tb)
		{
			index_t m = a.nrows();
			index_t k = a.ncolumns();
			index_t n = b.nrows();

			for (index_t j = 0; j < n; ++j)
			{
				for (index_t i = 0; i < m; ++i)
				{
					double v = 0;
					for (index_t l = 0; l < k; ++l) v += a(i, l) * b(j, l);
					r(i, j) = (T)(alpha * v + beta * c(i, j));
				}
			}
		}
		else if (ta && !tb)
		{
			index_t k = a.nrows();
			index_t m = a.ncolumns();
			index_t n = b.ncolumns();

			for (index_t j = 0; j < n; ++j)
			{
				for (index_t i = 0; i < m; ++i)
				{
					double v = 0;
					for (index_t l = 0; l < k; ++l) v += a(l, i) * b(l, j);
					r(i, j) = (T)(alpha * v + beta * c(i, j));
				}
			}
		}
		else // ta && tb
		{
			index_t k = a.nrows();
			index_t m = a.ncolumns();
			index_t n = b.nrows();

			for (index_t j = 0; j < n; ++j)
			{
				for (index_t i = 0; i < m; ++i)
				{
					double v = 0;
					for (index_t l = 0; l < k; ++l) v += a(l, i) * b(j, l);
					r(i, j) = (T)(alpha * v + beta * c(i, j));
				}
			}
		}
	}


} }

#endif

