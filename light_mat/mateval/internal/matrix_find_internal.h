/**
 * @file matrix_find_internal.h
 *
 * @brief Internal implementation of matrix-find algorithms
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_FIND_INTERNAL_H_
#define LIGHTMAT_MATRIX_FIND_INTERNAL_H_

#include <light_mat/mateval/ewise_eval.h>

namespace lmat { namespace internal {

	// count

	template<bool IsLinear> struct count_impl;

	template<class Reader>
	LMAT_ENSURE_INLINE
	inline size_t count_by_reader(index_t n, const Reader& rd)
	{
		size_t cnt = 0;
		for (index_t i = 0; i < n; ++i)
		{
			if (rd.scalar(i)) ++cnt;
		}
		return cnt;
	}

	template<>
	struct count_impl<true>
	{
		template<class A>
		LMAT_ENSURE_INLINE
		static size_t run(const A& a)
		{
			return count_by_reader(a.nelems(),
					make_vec_accessor(scalar_(), in_(a)));

		}
	};

	template<>
	struct count_impl<false>
	{
		template<class A>
		static size_t run(const A& a)
		{
			auto rd = make_multicol_accessor(scalar_(), in_(a));
			const index_t m = a.nrows();
			const index_t n = a.ncolumns();

			size_t cnt = 0;
			for (index_t j = 0; j < n; ++j)
				cnt += count_by_reader(m, rd.col(j));

			return cnt;
		}
	};

} }

#endif /* MATRIX_ALGS_INTERNAL_H_ */
