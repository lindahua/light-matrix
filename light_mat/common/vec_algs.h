/**
 * @file vec_algs.h
 *
 * @brief Simple algorithms on vectors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VEC_ALGS_H_
#define LIGHTMAT_VEC_ALGS_H_

#include <light_mat/common/basic_defs.h>
#include <functional>
#include <cmath>

namespace lmat
{

	// vec_all

	template<typename Pred, typename T>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const T *a)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i])) return false;
		}
		return true;
	}

	template<typename Pred, typename T>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const T *a, index_t astep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i * astep])) return false;
		}
		return true;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const Ta *a, const Tb *b)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i], b[i])) return false;
		}
		return true;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const Ta *a, index_t astep, const Tb *b)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i * astep], b[i])) return false;
		}
		return true;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const Ta *a, const Tb *b, index_t bstep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i], b[i * bstep])) return false;
		}
		return true;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_all(Pred pred, const index_t n, const Ta *a, index_t astep, const Tb *b, index_t bstep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (!pred(a[i * astep], b[i * bstep])) return false;
		}
		return true;
	}


	// vec_any

	template<typename Pred, typename T>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const T *a)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i])) return true;
		}
		return false;
	}

	template<typename Pred, typename T>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const T *a, index_t astep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i * astep])) return true;
		}
		return false;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const Ta *a, const Tb *b)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i], b[i])) return true;
		}
		return false;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const Ta *a, index_t astep, const Tb *b)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i * astep], b[i])) return true;
		}
		return false;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const Ta *a, const Tb *b, index_t bstep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i], b[i * bstep])) return true;
		}
		return false;
	}

	template<typename Pred, typename Ta, typename Tb>
	LMAT_ENSURE_INLINE
	bool vec_any(Pred pred, const index_t n, const Ta *a, index_t astep, const Tb *b, index_t bstep)
	{
		for (index_t i = 0; i < n; ++i)
		{
			if (pred(a[i * astep], b[i * bstep])) return true;
		}
		return false;
	}


	// vec_equal

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, const T& v)
	{
		return vec_all(std::bind2nd(std::equal_to<T>(), v), n, a);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, index_t astep, const T& v)
	{
		return vec_all(std::bind2nd(std::equal_to<T>(), v), n, a, astep);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, const T *b)
	{
		return vec_all(std::equal_to<T>(), n, a, b);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, index_t astep, const T *b)
	{
		return vec_all(std::equal_to<T>(), n, a, astep, b);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, const T *b, index_t bstep)
	{
		return vec_all(std::equal_to<T>(), n, a, b, bstep);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_equal(const index_t n, const T *a, index_t astep, const T *b, index_t bstep)
	{
		return vec_all(std::equal_to<T>(), n, a, astep, b, bstep);
	}


	// vec_approx

	LMAT_ENSURE_INLINE
	inline bool is_approx(double x, double y, double tol)
	{
		return std::fabs(x - y) < tol;
	}

	LMAT_ENSURE_INLINE
	inline bool is_approx(float x, float y, float tol)
	{
		return std::fabs(x - y) < tol;
	}

	template<typename T>
	struct approx_to : public std::binary_function<T, T, bool>
	{
		const T tolerance;

		LMAT_ENSURE_INLINE
		approx_to(const T& tol) : tolerance(tol) { }

		LMAT_ENSURE_INLINE
		bool operator() (const T& x, const T& y) const
		{
			return is_approx(x, y, tolerance);
		}
	};


	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, const T& v, const T& tol)
	{
		return vec_all(std::bind2nd(approx_to<T>(tol), v), n, a);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, index_t astep, const T& v, const T& tol)
	{
		return vec_all(std::bind2nd(approx_to<T>(tol), v), n, a, astep);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, const T *b, const T& tol)
	{
		return vec_all(approx_to<T>(tol), n, a, b);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, index_t astep, const T *b, const T& tol)
	{
		return vec_all(approx_to<T>(tol), n, a, astep, b);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, const T *b, index_t bstep, const T& tol)
	{
		return vec_all(approx_to<T>(tol), n, a, b, bstep);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	bool vec_approx(const index_t n, const T *a, index_t astep, const T *b, index_t bstep, const T& tol)
	{
		return vec_all(approx_to<T>(tol), n, a, astep, b, bstep);
	}

}

#endif /* VEC_ALGS_H_ */



