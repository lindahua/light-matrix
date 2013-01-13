/**
 * @file mat_fold.h
 *
 * Matrix folding (this is the basis for matrix reduction)
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_FOLD_H_
#define LIGHTMAT_MAT_FOLD_H_

#include <light_mat/mateval/ewise_eval.h>
#include "internal/mat_fold_internal.h"

#include <utility>


#define LMAT_DEFINE_SIMPLE_FOLDER( Name, InitExpr, FoldExpr, ReducExpr ) \
	template<typename T> \
	struct Name##_folder { \
		typedef T value_type; \
		LMAT_ENSURE_INLINE \
		T init(const T& x) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		void fold(T& a, const T& x) const { FoldExpr; } \
	}; \
	template<typename T, typename Kind> \
	struct Name##_folder<simd_pack<T, Kind> > { \
		typedef simd_pack<T, Kind> value_type; \
		LMAT_ENSURE_INLINE \
		value_type init(const value_type& x) const { return InitExpr; } \
		LMAT_ENSURE_INLINE \
		void fold(value_type& a, const value_type& x) const { FoldExpr; } \
		LMAT_ENSURE_INLINE \
		T reduce(const value_type& a) const { return ReducExpr; } \
	}; \
	LMAT_DECL_SIMDIZABLE_ON_REAL( Name##_folder ) \
	LMAT_DEF_TRIVIAL_SIMDIZE_MAP( Name##_folder )


namespace lmat
{

	/********************************************
	 *
	 *  common folders
	 *
	 ********************************************/

	LMAT_DEFINE_SIMPLE_FOLDER( sum, x, a += x, sum(a) )

	LMAT_DEFINE_SIMPLE_FOLDER( maximum, x, a = math::max(a, x), maximum(a) )

	LMAT_DEFINE_SIMPLE_FOLDER( minimum, x, a = math::min(a, x), minimum(a) )



	/********************************************
	 *
	 *  vectorized folding kernel
	 *
	 ********************************************/

	template<class Folder, typename U>
	class vecfold_kernel;

	template<class Folder, typename Kind>
	class vecfold_kernel<Folder, atags::simd<Kind> >
	{
		static_assert(is_simdizable<Folder, Kind>::value, "Folder should supports SIMD");

	public:
		typedef typename Folder::value_type value_type;

		LMAT_ENSURE_INLINE
		vecfold_kernel(const Folder& folder)
		: m_folder(folder) { }

		template<int N, typename Accessor>
		LMAT_ENSURE_INLINE
		value_type apply(dimension<N> dim, const Accessor& acc)
		{
			return internal::fold_impl(dim, atags::simd<Kind>(), m_folder, acc);
		}

		template<typename Accessor>
		LMAT_ENSURE_INLINE
		value_type apply(index_t len, const Accessor& acc)
		{
			return apply(dimension<0>(len), acc);
		}

		template<int N, typename Wrap>
		LMAT_ENSURE_INLINE
		value_type operator() (dimension<N> dim, const Wrap& wrap)
		{
			return apply(dim, make_vec_accessor(atags::simd<Kind>(), wrap));
		}

		template<int N, typename Wrap>
		LMAT_ENSURE_INLINE
		value_type operator() (index_t len, const Wrap& wrap)
		{
			return apply(dimension<0>(len), make_vec_accessor(atags::simd<Kind>(), wrap));
		}

	private:
		Folder m_folder;
	};


	/********************************************
	 *
	 *  convenient function
	 *
	 ********************************************/

	template<class Folder, typename U>
	LMAT_ENSURE_INLINE
	inline vecfold_kernel<Folder, U> fold(const Folder& folder, U)
	{
		return vecfold_kernel<Folder, U>(folder);
	}

	template<class Folder>
	LMAT_ENSURE_INLINE
	inline vecfold_kernel<Folder, default_access_unit_t> fold(const Folder& folder)
	{
		return vecfold_kernel<Folder, default_access_unit_t>(folder);
	}


}

#endif 
