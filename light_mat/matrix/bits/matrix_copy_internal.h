/*
 * @file matrix_copy_internal.h
 *
 * Internal implementation of matrix copy
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_COPY_INTERNAL_H_
#define LIGHTMAT_MATRIX_COPY_INTERNAL_H_

#include <light_mat/core/mem_op.h>

namespace lmat { namespace detail {

	template<typename T, int M, int N>
	struct fixed_mat_copy
	{
		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, T *dst)
		{
			copy_mem(M * N, src, dst);
		}


		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, T *dst, const index_t ldim_d)
		{
			for (index_t j = 0; j < N; ++j, src+=M, dst+=ldim_d)
			{
				copy_mem(M, src, dst);
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, const index_t ldim_s, T *dst)
		{
			for (index_t j = 0; j < N; ++j, src+=ldim_s, dst+=M)
			{
				copy_mem(M, src, dst);
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, const index_t ldim_s,
				      T *dst, const index_t ldim_d)
		{
			for (index_t j = 0; j < N; ++j, src+=ldim_s, dst+=ldim_d)
			{
				copy_mem(M, src, dst);
			}
		}
	};


	template<typename T>
	struct scalar_copy
	{
		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, T *dst)
		{
			*dst = *src;
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, T *dst, const index_t)
		{
			*dst = *src;
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, const index_t, T *dst)
		{
			*dst = *src;
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t,
				const T *src, const index_t,
				      T *dst, const index_t)
		{
			*dst = *src;
		}
	};


	template<typename T>
	struct col_copy
	{
		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t,
				const T *src, T *dst)
		{
			copy_mem(m, src, dst);
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t,
				const T *src, T *dst, const index_t)
		{
			copy_mem(m, src, dst);
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t,
				const T *src, const index_t, T *dst)
		{
			copy_mem(m, src, dst);
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t,
				const T *src, const index_t,
				      T *dst, const index_t)
		{
			copy_mem(m, src, dst);
		}
	};

	template<typename T>
	struct row_copy
	{
		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t n,
				const T *src, T *dst)
		{
			copy_mem(n, src, dst);
		}


		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t n,
				const T *src, T *dst, const index_t ldim_d)
		{
			for (index_t j = 0; j < n; ++j)
			{
				dst[j * ldim_d] = src[j];
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t n,
				const T *src, const index_t ldim_s, T *dst)
		{
			for (index_t j = 0; j < n; ++j)
			{
				dst[j] = src[j * ldim_s];
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t, const index_t n,
				const T *src, const index_t ldim_s,
				      T *dst, const index_t ldim_d)
		{
			for (index_t j = 0; j < n; ++j)
			{
				dst[j * ldim_d] = src[j * ldim_s];
			}
		}
	};

	template<typename T>
	struct mat_copy
	{
		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t n,
				const T *src, T *dst)
		{
			copy_mem(m * n, src, dst);
		}


		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t n,
				const T *src, T *dst, const index_t ldim_d)
		{
			if (n == 1)
			{
				copy_mem(m, src, dst);
			}
			else if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					dst[j * ldim_d] = src[j];
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j,
					src += m, dst += ldim_d)
				{
					copy_mem(m, src, dst);
				}
			}
		}

		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t n,
				const T *src, const index_t ldim_s, T *dst)
		{
			if (n == 1)
			{
				copy_mem(m, src, dst);
			}
			else if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					dst[j] = src[j * ldim_s];
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j,
					src += ldim_s, dst += m)
				{
					copy_mem(m, src, dst);
				}
			}
		}


		LMAT_ENSURE_INLINE
		static void copy(const index_t m, const index_t n,
				const T *src, const index_t ldim_s,
					  T *dst, const index_t ldim_d)
		{
			if (n == 1)
			{
				copy_mem(m, src, dst);
			}
			else if (m == 1)
			{
				for (index_t j = 0; j < n; ++j)
				{
					dst[j * ldim_d] = src[j * ldim_s];
				}
			}
			else
			{
				for (index_t j = 0; j < n; ++j,
					src += ldim_s, dst += ldim_d)
				{
					copy_mem(m, src, dst);
				}
			}
		}
	};


	template<typename T, int M, int N>
	struct mat_copier
	{
		typedef typename
				if_c<(N == 1),
					typename
					if_c<(M == 1),
						scalar_copy<T>, // N == 1 && M == 1
						col_copy<T> 	// N == 1 && M != 1
					>::type,
					typename
					if_c<(M == 1),
						row_copy<T>,	// N != 1 && M == 1
						typename
						if_c<((M > 0) && (N > 0)),
							fixed_mat_copy<T, M, N>, 	// M > 1 && N > 1
							mat_copy<T> 				// otherwise
						>::type
					>::type
				>::type type;
	};

} }

#endif 
