/**
 * @file mat_allany.h
 *
 * @brief all and any reduction on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_ALLANY_H_
#define LIGHTMAT_MAT_ALLANY_H_

#include "internal/mat_allany_internal.h"
#include <light_mat/mateval/mat_pred.h>

namespace lmat
{

	/********************************************
	 *
	 *  full reduction
	 *
	 ********************************************/

	template<typename T, class Mat>
	inline bool all(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		static_assert(internal::prefers_linear<Mat>::value,
				"mat should allow linear accessing");

		typedef default_simd_kind kind;
		const bool use_simd = internal::prefers_simd<Mat, T, kind, true>::value;

		typedef typename meta::if_c<use_simd, atags::simd<kind>, atags::scalar>::type atag;

		dimension<meta::nelems<Mat>::value> dim(mat.nelems());

		if (val)
		{
			return internal::all_impl(dim, type_<T>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
		else
		{
			return !internal::any_impl(dim, type_<T>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
	}


	template<class Mat>
	inline bool all(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		static_assert(internal::prefers_linear<Mat>::value,
				"mat should allow linear accessing");

		typedef atags::scalar atag;
		dimension<meta::nelems<Mat>::value> dim(mat.nelems());

		if (val)
		{
			return internal::all_impl(dim, type_<bool>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
		else
		{
			return !internal::any_impl(dim, type_<bool>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
	}


	template<typename T, class Mat>
	inline bool any(const IEWiseMatrix<Mat, mask_t<T> >& mat, bool val=true)
	{
		static_assert(internal::prefers_linear<Mat>::value,
				"mat should allow linear accessing");

		typedef default_simd_kind kind;
		const bool use_simd = internal::prefers_simd<Mat, T, kind, true>::value;

		typedef typename meta::if_c<use_simd, atags::simd<kind>, atags::scalar>::type atag;

		dimension<meta::nelems<Mat>::value> dim(mat.nelems());

		if (val)
		{
			return internal::any_impl(dim, type_<T>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
		else
		{
			return !internal::all_impl(dim, type_<T>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
	}

	template<class Mat>
	inline bool any(const IEWiseMatrix<Mat, bool>& mat, bool val=true)
	{
		static_assert(internal::prefers_linear<Mat>::value,
				"mat should allow linear accessing");

		typedef atags::scalar atag;
		dimension<meta::nelems<Mat>::value> dim(mat.nelems());

		if (val)
		{
			return internal::any_impl(dim, type_<bool>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
		else
		{
			return !internal::all_impl(dim, type_<bool>(), atag(),
					make_vec_accessor(atag(), in_(mat.derived())));
		}
	}

}

#endif
