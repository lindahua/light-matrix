/**
 * @file matrix_indexer.h
 *
 * Matrix entry index calculation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_INDEXER_H_
#define LIGHTMAT_MATRIX_INDEXER_H_

#include <light_mat/core/basic_defs.h>

namespace lmat
{
	namespace detail
	{
		template<bool IsRow, bool IsCol> struct matrix_index_helper;

		template<> struct matrix_index_helper<false, false>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return i + j * ldim;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return i + j * ldim;
			}
		};

		template<> struct matrix_index_helper<false, true>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return i;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return i;
			}
		};

		template<> struct matrix_index_helper<true, false>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return j * ldim;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return j;
			}
		};

		template<> struct matrix_index_helper<true, true>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return 0;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return 0;
			}
		};
	}


	template<int M, int N>
	struct matrix_indexer
	{
		LMAT_ENSURE_INLINE
		static index_t offset(const index_t ldim, const index_t i, const index_t j)
		{
			return detail::matrix_index_helper<M == 1, N == 1>::offset(ldim, i, j);
		}

		LMAT_ENSURE_INLINE
		static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
		{
			return detail::matrix_index_helper<M == 1, N == 1>::offset_c(ldim, i, j);
		}
	};


}

#endif /* OFFSET_CALC_H_ */
