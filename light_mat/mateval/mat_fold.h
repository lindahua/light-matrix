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
		void fold(T& a, const T& b) const { a += b; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& b) const { a += b; }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& p) const { return math::sum(p); }
	};

	template<typename T>
	struct folder_supports_simd<sum_folder<T> >
	{
		static const bool value = true;
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
		void fold(T& a, const T& b) const { a = math::max(a, b); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& b) const { a = math::max(a, b); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& p) const { return math::maximum(p); }
	};

	template<typename T>
	struct folder_supports_simd<maximum_folder<T> >
	{
		static const bool value = true;
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
		void fold(T& a, const T& b) const { a = math::min(a, b); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void fold(math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& b) const { a = math::min(a, b); }

		template<typename Kind>
		LMAT_ENSURE_INLINE
		T reduce(const math::simd_pack<T, Kind>& p) const { return math::minimum(p); }
	};

	template<typename T>
	struct folder_supports_simd<minimum_folder<T> >
	{
		static const bool value = true;
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
		typedef typename meta::fun_value_type<Folder>::type value_type;

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


	template<class Folder, class TermFun, typename U>
	class vecfoldf_kernel
	{
	public:
		typedef typename meta::fun_value_type<Folder>::type value_type;

		LMAT_ENSURE_INLINE
		vecfoldf_kernel(const Folder& folder, const TermFun& termf)
		: m_folder(folder), m_termfun(termf) { }

		template<int N, typename... Accessors>
		LMAT_ENSURE_INLINE
		value_type apply(dimension<N> dim, const Accessors&... accs)
		{
			return internal::foldf_impl(dim, U(), m_folder, m_termfun, accs...);
		}

		template<typename... Accessors>
		LMAT_ENSURE_INLINE
		value_type apply(index_t len, const Accessors&... accs)
		{
			return apply(dimension<0>(len), accs...);
		}

		template<int N, typename... Wraps>
		LMAT_ENSURE_INLINE
		value_type operator() (dimension<N> dim, const Wraps&... wraps)
		{
			return apply(dim, make_vec_accessor(U(), wraps)...);
		}

		template<int N, typename... Wraps>
		LMAT_ENSURE_INLINE
		value_type operator() (index_t len, const Wraps&... wraps)
		{
			return apply(dimension<0>(len), make_vec_accessor(U(), wraps)...);
		}

	private:
		Folder m_folder;
		TermFun m_termfun;
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

	template<class Folder, class TermFun, typename U>
	LMAT_ENSURE_INLINE
	inline vecfoldf_kernel<Folder, TermFun, U> foldf(const Folder& folder, const TermFun& tfun, U)
	{
		return vecfoldf_kernel<Folder, TermFun, U>(folder, tfun);
	}

	template<class Folder, class TermFun>
	LMAT_ENSURE_INLINE
	inline vecfoldf_kernel<Folder, TermFun, default_access_unit_t> foldf(const Folder& folder, const TermFun& tfun)
	{
		return vecfoldf_kernel<Folder, TermFun, default_access_unit_t>(folder, tfun);
	}


}

#endif 
