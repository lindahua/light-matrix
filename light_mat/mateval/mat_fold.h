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

	/********************************************
	 *
	 *  useful folders
	 *
	 ********************************************/

	template<typename T>
	struct sum_folder
	{
		typedef T value_type;

		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::sum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a += b;
		}
	};

	template<typename T>
	struct maximum_folder
	{
		typedef T value_type;

		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::maximum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a = math::max(a, b);
		}
	};

	template<typename T>
	struct minimum_folder
	{
		typedef T value_type;

		template<class Pack>
		LMAT_ENSURE_INLINE
		T reduce(const Pack& p) const
		{
			return math::minimum(p);
		}

		template<class A>
		LMAT_ENSURE_INLINE
		void fold(A& a, const A& b) const
		{
			a = math::min(a, b);
		}
	};


	/********************************************
	 *
	 *  vectorized folding kernel
	 *
	 ********************************************/

	template<class Folder, typename U>
	class vecfold_kernel
	{
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
