/**
 * @file matrix_ordstats.h
 *
 * Algorithms on matrices related to order statistics
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ORDSTATS_H_
#define LIGHTMAT_MATRIX_ORDSTATS_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/ewise_eval.h>
#include <utility>
#include <algorithm>

namespace lmat
{
	/********************************************
	 *
	 *  finding max/min
	 *
	 ********************************************/

	namespace internal
	{
		template<class A, typename T>
		inline void _find_max(const IEWiseMatrix<A, T>& a, index_t& p, T& s)
		{
			const index_t n = a.nelems();

			auto rd = make_vec_accessor(atags::scalar(), in_(a.derived()));
			p = 0;
			s = rd.scalar(0);

			for (index_t i = 1; i < n; ++i)
			{
				T x = rd.scalar(i);
				if (x > a)
				{
					p = i;
					s = x;
				}
			}

			return p;
		}

		template<class A, typename T>
		inline void _find_min(const IEWiseMatrix<A, T>& a, index_t& p, T& s)
		{
			const index_t n = a.nelems();

			auto rd = make_vec_accessor(atags::scalar(), in_(a.derived()));
			p = 0;
			s = rd.scalar(0);

			for (index_t i = 1; i < n; ++i)
			{
				T x = rd.scalar(i);
				if (x < a)
				{
					p = i;
					s = x;
				}
			}

			return p;
		}
	}

	template<class A, typename T>
	inline std::enable_if<supports_linear_macc<A>::value, index_t>::type
	find_imax(const IEWiseMatrix<A, T>& a)
	{
		if (is_empty(a))
			throw invalid_argument("find_imax: argument a was empty.");

		index_t p;
		T s;
		internal::_find_max(a, p, s);
		return p;
	}

	template<class A, typename T>
	inline std::enable_if<supports_linear_macc<A>::value, index_t>::type
	find_imin(const IEWiseMatrix<A, T>& a)
	{
		if (is_empty(a))
			throw invalid_argument("find_imin: argument a was empty.");

		index_t p;
		T s;
		internal::_find_min(a, p, s);
		return p;
	}

	template<class A, typename T>
	inline std::enable_if<supports_linear_macc<A>::value, std::pair<index_t, T> >::type
	find_max(const IEWiseMatrix<A, T>& a)
	{
		if (is_empty(a))
			throw invalid_argument("find_max: argument a was empty.");

		index_t p;
		T s;
		internal::_find_max(a, p, s);
		return std::make_pair(p, s);
	}

	template<class A, typename T>
	inline std::enable_if<supports_linear_macc<A>::value, std::pair<index_t, T> >::type
	find_min(const IEWiseMatrix<A, T>& a)
	{
		if (is_empty(a))
			throw invalid_argument("find_min: argument a was empty.");

		index_t p;
		T s;
		internal::_find_min(a, p, s);
		return std::make_pair(p, s);
	}


	template<typename T, class A, typename TI, class D>
	inline std::enable_if<meta::supports_linear_index<D>::value, void>::type
	colwise_find_imax(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TI>& idx)
	{
		if (is_empty(a))
			throw invalid_argument("colwise_find_imax: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();

		const index_t n = a.ncolumns();
		index_t p;
		T s;

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_max(a.column(j), p, s);
			idx_[j] = static_cast<TI>(p);
		}
	}

	template<typename T, class A, typename TI, class D>
	inline typename std::enable_if<meta::supports_linear_index<D>::value,
	void>::type
	colwise_find_imin(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TI>& idx)
	{
		if (is_empty(a))
			throw invalid_argument("colwise_find_imin: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();

		const index_t n = a.ncolumns();
		index_t p;
		T s;

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_min(a.column(j), p, s);
			idx_[j] = static_cast<TI>(p);
		}
	}


	template<typename T, class A, typename TI, class D, class R>
	inline typename std::enable_if<
		meta::supports_linear_index<D>::value &&
		meta::supports_linear_index<R>::value,
	void>::type
	colwise_find_max(const IEWiseMatrix<A, T>& a,
			IRegularMatrix<D, TI>& idx, IRegularMatrix<R, T>& r)
	{
		if (is_empty(a))
			throw invalid_argument("colwise_find_max: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();
		R& r_ = r.derived();

		const index_t n = a.ncolumns();
		index_t p;
		T s;

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_max(a.column(j), p, s);
			idx_[j] = static_cast<TI>(p);
			r_[j] = s;
		}
	}

	template<typename T, class A, typename TI, class D, class R>
	inline typename std::enable_if<
		meta::supports_linear_index<D>::value &&
		meta::supports_linear_index<R>::value,
	void>::type
	colwise_find_min(const IEWiseMatrix<A, T>& a,
			IRegularMatrix<D, TI>& idx, IRegularMatrix<R, T>& r)
	{
		if (is_empty(a))
			throw invalid_argument("colwise_find_min: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();
		R& r_ = r.derived();

		const index_t n = a.ncolumns();
		index_t p;
		T s;

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_min(a.column(j), p, s);
			idx_[j] = static_cast<TI>(p);
			r_[j] = s;
		}
	}


	/********************************************
	 *
	 *  finding elements at specific location
	 *
	 ********************************************/

	namespace internal
	{
		template<typename T, class A>
		T _find_nth(IRegularMatrix<A, T>& a, index_t k)
		{
			T *p = a.ptr_data();
			index_t n = a.elems();

			std::nth_element(p, p+k, p+n);
			return p[k];
		}
	}


}

#endif 
