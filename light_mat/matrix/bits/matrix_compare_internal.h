/**
 * @file matrix_compare_internal.h
 *
 * Internal implementation of matrix comparison
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COMPARE_INTERNAL_H_
#define LIGHTMAT_MATRIX_COMPARE_INTERNAL_H_

#include <light_mat/core/mem_op.h>

namespace lmat { namespace detail {

	template<typename T, int M, int N>
	struct fixed_mat_compare
	{
		LMAT_ENSURE_INLINE
		static bool is_equal(const index_t, const index_t,
				const T *a, const index_t ldim_a,
				const T *b, const index_t ldim_b)
		{
			if (ldim_a == M && ldim_b == M)
			{
				return mem_equal(M * N, a, b);
			}
			else
			{
				for (index_t j = 0; j < N; ++j, a+=ldim_a, b+=ldim_b)
				{
					if (!mem_equal(M, a, b)) return false;
				}
				return true;
			}
		}
	};

	template<typename T>
	struct scalar_compare
	{
		LMAT_ENSURE_INLINE
		static bool is_equal(const index_t, const index_t,
				const T *a, const index_t,
				const T *b, const index_t)
		{
			return *a == *b;
		}
	};

	template<typename T>
	struct col_compare
	{
		LMAT_ENSURE_INLINE
		static bool is_equal(const index_t m, const index_t,
				const T *a, const index_t,
				const T *b, const index_t)
		{
			return mem_equal(m, a, b);
		}
	};

	template<typename T>
	struct row_compare
	{
		LMAT_ENSURE_INLINE
		static bool is_equal(const index_t, const index_t n,
				const T *a, const index_t ldim_a,
				const T *b, const index_t ldim_b)
		{
			if (ldim_a == 1 && ldim_b == 1)
			{
				return mem_equal(n, a, b);
			}
			else
			{
				for (index_t j = 0; j < n; ++j)
				{
					if (!(a[j * ldim_a] == b[j * ldim_b])) return false;
				}
				return true;
			}
		}
	};

	template<typename T>
	struct mat_compare
	{
		inline
		static bool is_equal(const index_t m, const index_t n,
				const T *a, const index_t ldim_a,
				const T *b, const index_t ldim_b)
		{
			bool eq;

			if (n == 1)
			{
				eq = mem_equal(m, a, b);
			}
			else if (ldim_a == m && ldim_b == m)
			{
				eq = mem_equal(m * n, a, b);
			}
			else if (m == 1)
			{
				eq = true;
				for (index_t j = 0; j < n; ++j)
				{
					if ( !(a[j * ldim_a] == b[j * ldim_b]) )
					{
						eq = false;
						break;
					}
				}
			}
			else
			{
				eq = true;
				for (index_t j = 0; j < n; ++j, a += ldim_a, b += ldim_b)
				{
					if ( !mem_equal(m, a, b) )
					{
						eq = false;
						break;
					}
				}
			}

			return eq;
		}
	};


	template<typename T, int M, int N>
	struct mat_comparer
	{
		typedef typename
				if_c<(N == 1),
					typename
					if_c<(M == 1),
						scalar_compare<T>, // N == 1 && M == 1
						col_compare<T> 	// N == 1 && M != 1
					>::type,
					typename
					if_c<(M == 1),
						row_compare<T>,	// N != 1 && M == 1
						typename
						if_c<((M > 0) && (N > 0)),
							fixed_mat_compare<T, M, N>, 	// M > 1 && N > 1
							mat_compare<T> 				// otherwise
						>::type
					>::type
				>::type type;
	};

} }

#endif /* MATRIX_COMPARE_INTERNAL_H_ */
