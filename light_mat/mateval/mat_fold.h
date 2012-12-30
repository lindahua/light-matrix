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

#ifndef LIGHTMAT_MAT_FOLD_H_
#define LIGHTMAT_MAT_FOLD_H_

#include <light_mat/mateval/ewise_eval.h>
#include "internal/mat_fold_internal.h"

#include <utility>

namespace lmat
{
	// folder interface

	template<class Folder>
	struct folder_supports_simd
	{
		static const bool value = false;
	};

	template<class Folder, typename Kind>
	struct folder_simd_pack;


	/********************************************
	 *
	 *  common folders
	 *
	 ********************************************/

	// sum

	template<typename T>
	struct sum_folder
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		T init(const T& x) const { return x; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		math::simd_pack<T, Kind> init(const math::simd_pack<T, Kind>& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		void fold(T& a, const T& x) const { a += x; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& x) const { a += x; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& a) const { return math::sum(a); }
	};

	template<typename T>
	struct folder_supports_simd<sum_folder<T> >
	{
		static const bool value =
				std::is_same<T, float>::value || std::is_same<T, double>::value;
	};

	template<typename T, typename Kind>
	struct folder_simd_pack<sum_folder<T>, Kind>
	{
		typedef math::simd_pack<T, Kind> type;
	};

	// maximum

	template<typename T>
	struct maximum_folder
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		T init(const T& x) const { return x; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		math::simd_pack<T, Kind> init(const math::simd_pack<T, Kind>& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		void fold(T& a, const T& x) const { a = math::max(a, x); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& x) const { a = math::max(a, x); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& a) const { return math::maximum(a); }
	};

	template<typename T>
	struct folder_supports_simd<maximum_folder<T> >
	{
		static const bool value =
				std::is_same<T, float>::value || std::is_same<T, double>::value;
	};

	template<typename T, typename Kind>
	struct folder_simd_pack<maximum_folder<T>, Kind>
	{
		typedef math::simd_pack<T, Kind> type;
	};

	// minimum

	template<typename T>
	struct minimum_folder
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		T init(const T& x) const { return x; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		math::simd_pack<T, Kind> init(const math::simd_pack<T, Kind>& x) const
		{
			return x;
		}

		LMAT_ENSURE_INLINE
		void fold(T& a, const T& x) const { a = math::min(a, x); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& x) const { a = math::min(a, x); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& a) const { return math::minimum(a); }
	};

	template<typename T>
	struct folder_supports_simd<minimum_folder<T> >
	{
		static const bool value =
				std::is_same<T, float>::value || std::is_same<T, double>::value;
	};

	template<typename T, typename Kind>
	struct folder_simd_pack<minimum_folder<T>, Kind>
	{
		typedef math::simd_pack<T, Kind> type;
	};


	/********************************************
	 *
	 *  vectorized folding kernel
	 *
	 ********************************************/

	template<class Folder, typename U>
	class vecfold_kernel
	{
		static_assert(folder_supports_simd<Folder>::value, "Folder should supports SIMD");

	public:
		typedef typename Folder::value_type value_type;

		LMAT_ENSURE_INLINE
		vecfold_kernel(const Folder& folder)
		: m_folder(folder) { }

		template<int N, typename Accessor>
		LMAT_ENSURE_INLINE
		value_type apply(dimension<N> dim, const Accessor& acc)
		{
			return internal::fold_impl(dim, U(), m_folder, acc);
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
			return apply(dim, make_vec_accessor(U(), wrap));
		}

		template<int N, typename Wrap>
		LMAT_ENSURE_INLINE
		value_type operator() (index_t len, const Wrap& wrap)
		{
			return apply(dimension<0>(len), make_vec_accessor(U(), wrap));
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
