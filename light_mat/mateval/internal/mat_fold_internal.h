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

#include <light_mat/mateval/macc_policy.h>


namespace lmat { namespace internal {

	template<class FoldKernel, class Shape, typename... Args>
	struct fold_policy
	{
		static const bool supp_linear =
				meta::all_<supports_linear_access<Args>...>::value;

		typedef default_simd_kind skind;

		static const bool supp_simd =
				is_simdizable<FoldKernel, skind>::value &&
				meta::all_<supports_simd<Args, skind>...>::value;

		static const bool use_linear = supp_linear;

		static const index_t _len =
				use_linear ? (Shape::ct_nrows * Shape::ct_ncols) : Shape::ct_nrows;

		static const unsigned int pack_width =
				internal::_kernel_packwidth<FoldKernel, skind,
				is_simdizable<FoldKernel, skind>::value>::value;

		typedef typename std::conditional<use_linear,
				linear_, percol_>::type access;

		static const bool use_simd = supp_simd && ((unsigned int)_len % pack_width == 0);

		typedef typename std::conditional<use_simd, simd_<skind>, scalar_>::type unit;

		typedef macc_<access, unit> type;
	};

	/********************************************
	 *
	 *  core implementation
	 *
	 ********************************************/

	template<index_t Len, class FoldKernel, typename... Reader>
	LMAT_ENSURE_INLINE
	inline typename FoldKernel::accumulated_type
	linear_fold_impl(const dimension<Len>& dim, scalar_, const FoldKernel& fker, const Reader&... rd)
	{
		typedef typename FoldKernel::accumulated_type RT;
		RT r = fker.init(rd.scalar(0)...);
		const index_t len = dim.value();
		for (index_t i = 1; i < len; ++i) fker(r, rd.scalar(i)...);
		return r;
	}

	template<index_t Len, typename SKind, class FoldKernel, typename... Reader>
	inline typename FoldKernel::accumulated_type
	linear_fold_impl(const dimension<Len>& dim, simd_<SKind>, const FoldKernel& fker, const Reader&... rd)
	{
		typedef typename FoldKernel::accumulated_type RT;
		typedef typename simdize_map<FoldKernel, SKind>::type simd_fker_t;
		typedef typename simd_fker_t::accumulated_type pack_t;

		simd_fker_t pk_fker = simdize_map<FoldKernel, SKind>::get(fker);

		const index_t pw = (index_t)pack_t::pack_width;

		const index_t len = dim.value();
		const index_t npacks = (index_t)int_div<pack_t::pack_width>::quo((size_t)len);
		index_t i;
		RT r;

		if (npacks)
		{
			pass(rd.begin_packs()...);

			pack_t a0, a1, a2, a3;

			const index_t m4 = npacks >> 2;
			if (m4)
			{
				const index_t w4 = pw << 2;
				const index_t l4 = m4 * w4;

				a0 = pk_fker.init(rd.pack(0)...);
				a1 = pk_fker.init(rd.pack(pw)...);
				a2 = pk_fker.init(rd.pack(pw * 2)...);
				a3 = pk_fker.init(rd.pack(pw * 3)...);

				for (i = w4; i < l4; i += w4)
				{
					pk_fker(a0, rd.pack(i)...);
					pk_fker(a1, rd.pack(i + pw)...);
					pk_fker(a2, rd.pack(i + pw * 2)...);
					pk_fker(a3, rd.pack(i + pw * 3)...);
				}

				pk_fker(a0, a2);
				pk_fker(a1, a3);

				if (npacks & 2)
				{
					pk_fker(a0, rd.pack(i)...);
					pk_fker(a1, rd.pack(i + pw)...);
					i += pw * 2;
				}

				pk_fker(a0, a1);

				if (npacks & 1)
				{
					pk_fker(a0, rd.pack(i)...);
					i += pw;
				}
			}
			else if (npacks >> 1)
			{
				a0 = pk_fker.init(rd.pack(0)...);
				pk_fker(a0, rd.pack(pw)...);
				i = pw << 1;

				if (npacks & 1)
				{
					pk_fker(a0, rd.pack(i)...);
					i += pw;
				}
			}
			else
			{
				a0 = pk_fker.init(rd.pack(0)...);
				i = pw;
			}

			pass(rd.end_packs()...);

			r = pk_fker.reduce(a0);
		}
		else
		{
			r = fker.init(rd.scalar(0)...);
			i = 1;
		}

		for (; i < len; ++i) fker(r, rd.scalar(i)...);
		return r;
	}


	template<index_t CM, index_t CN, typename U, class FoldKernel, typename... Reader>
	inline typename FoldKernel::accumulated_type
	percol_fold_impl(const matrix_shape<CM, CN>& shape, U, const FoldKernel& fker, const Reader&... rd)
	{
		typedef typename FoldKernel::accumulated_type RT;

		dimension<CM> col_dim(shape.nrows());

		RT r = linear_fold_impl(col_dim, U(), fker, rd.col(0)...);

		const index_t n = shape.ncolumns();
		for (index_t j = 1; j < n; ++j)
		{
			RT rj = linear_fold_impl(col_dim, U(), fker, rd.col(j)...);
			fker(r, rj);
		}

		return r;
	}

} }

#endif 
