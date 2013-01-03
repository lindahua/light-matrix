/**
 * @file mat_fold_internal.h
 *
 * Internal implementation of vector folding
 * 
 * @author Dahua Lin 
 */


#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_FOLD_INTERNAL_H_
#define LIGHTMAT_MAT_FOLD_INTERNAL_H_

#include <light_mat/mateval/mateval_fwd.h>

namespace lmat {

	template<class Folder, typename Kind>
	struct folder_simd_pack;

	namespace internal {

	/********************************************
	 *
	 *  core implementation
	 *
	 ********************************************/

	template<int N, class Folder, class Reader>
	LMAT_ENSURE_INLINE
	inline typename Folder::value_type
	fold_impl(const dimension<N>& dim, atags::scalar, const Folder& folder, const Reader& rd)
	{
		typedef typename Folder::value_type T;
		T r = folder.init(rd.scalar(0));
		const index_t len = dim.value();
		for (index_t i = 1; i < len; ++i) folder.fold(r, rd.scalar(i));
		return r;
	}

	template<int N, typename SKind, class Folder, class Reader>
	inline typename Folder::value_type
	fold_impl(const dimension<N>& dim, atags::simd<SKind>, const Folder& folder, const Reader& rd)
	{
		typedef typename Folder::value_type T;
		typedef typename folder_simd_pack<Folder, SKind>::type pack_t;

		const index_t pw = (index_t)pack_t::pack_width;

		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i;
		T r;

		if (npacks)
		{
			rd.begin_packs();

			pack_t a0, a1, a2, a3;

			const index_t m4 = npacks >> 2;
			if (m4)
			{
				const index_t w4 = pw << 2;
				const index_t l4 = m4 * w4;

				a0 = folder.init(rd.pack(0));
				a1 = folder.init(rd.pack(pw));
				a2 = folder.init(rd.pack(pw * 2));
				a3 = folder.init(rd.pack(pw * 3));

				for (i = w4; i < l4; i += w4)
				{
					folder.fold(a0, rd.pack(i));
					folder.fold(a1, rd.pack(i + pw));
					folder.fold(a2, rd.pack(i + pw * 2));
					folder.fold(a3, rd.pack(i + pw * 3));
				}

				folder.fold(a0, a2);
				folder.fold(a1, a3);

				if (npacks & 2)
				{
					folder.fold(a0, rd.pack(i));
					folder.fold(a1, rd.pack(i + pw));
					i += pw * 2;
				}

				folder.fold(a0, a1);

				if (npacks & 1)
				{
					folder.fold(a0, rd.pack(i));
					i += pw;
				}
			}
			else if (npacks >> 1)
			{
				a0 = folder.init(rd.pack(0));
				folder.fold(a0, rd.pack(pw));
				i = pw << 1;

				if (npacks & 1)
				{
					folder.fold(a0, rd.pack(i));
					i += pw;
				}
			}
			else
			{
				a0 = folder.init(rd.pack(0));
				i = pw;
			}

			rd.end_packs();

			r = folder.reduce(a0);
		}
		else
		{
			r = folder.init(rd.scalar(0));
			i = 1;
		}

		for (; i < len; ++i) folder.fold(r, rd.scalar(i));
		return r;
	}

} }

#endif 
