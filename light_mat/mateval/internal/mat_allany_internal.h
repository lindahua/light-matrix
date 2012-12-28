/**
 * @file mat_allany_internal.h
 *
 * @brief Internal implementation of all/any reduction
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_ALLANY_INTERNAL_H_
#define LIGHTMAT_MAT_ALLANY_INTERNAL_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/math/sse_reduce.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_reduce.h>
#endif

namespace lmat { namespace internal {


	template<int N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool all_impl(const dimension<N>& dim, type_<T>, atags::scalar, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (!rd.scalar(i)) return false;
		}

		return true;
	}


	template<int N, typename T, typename SKind, class Reader>
	inline bool all_impl(const dimension<N>& dim, type_<T>, atags::simd<SKind>, const Reader& rd)
	{
		typedef typename math::simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = int_div<pack_t::pack_width>::quo(len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (math::any_false(rd.pack(i))) return false;
			}
		}

		return all_impl(dim, type_<T>(), atags::scalar(), rd, i);
	}


	template<int N, typename T, class Reader>
	LMAT_ENSURE_INLINE
	inline bool any_impl(const dimension<N>& dim, type_<T>, atags::scalar, const Reader& rd, index_t i=0)
	{
		const index_t len = dim.value();

		for (; i < len; ++i)
		{
			if (rd.scalar(i)) return true;
		}

		return false;
	}


	template<int N, typename T, typename SKind, class Reader>
	inline bool any_impl(const dimension<N>& dim, type_<T>, atags::simd<SKind>, const Reader& rd)
	{
		typedef typename math::simd_bpack<T, SKind> pack_t;

		const index_t pw = (index_t)pack_t::pack_width;
		const index_t len = dim.value();
		const index_t npacks = int_div<pack_t::pack_width>::quo(len);
		index_t i = 0;

		if (npacks)
		{
			const index_t l = npacks * pw;
			for (i = 0; i < l; i += pw)
			{
				if (math::any_true(rd.pack(i))) return true;
			}
		}

		return any_impl(dim, type_<T>(), atags::scalar(), rd, i);
	}


} }

#endif /* MAT_ALLANY_INTERNAL_H_ */
