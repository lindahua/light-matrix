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
		template<class Rd, typename T>
		inline void _find_max(const Rd& rd, index_t n, index_t& p, T& s)
		{
			p = 0;
			s = rd.scalar(0);

			for (index_t i = 1; i < n; ++i)
			{
				T x = rd.scalar(i);
				if (x > s)
				{
					p = i;
					s = x;
				}
			}
		}

		template<class Rd, typename T>
		inline void _find_min(const Rd& rd, index_t n, index_t& p, T& s)
		{
			p = 0;
			s = rd.scalar(0);

			for (index_t i = 1; i < n; ++i)
			{
				T x = rd.scalar(i);
				if (x < s)
				{
					p = i;
					s = x;
				}
			}
		}
	}

	template<class A, typename T>
	inline typename std::enable_if<supports_linear_macc<A>::value,
	index_t>::type
	find_imax(const IEWiseMatrix<A, T>& a)
	{
		const index_t n = a.nelems();
		if (n == 0)
			throw invalid_argument("find_imax: argument a was empty.");

		index_t p;
		T s;
		internal::_find_max(make_vec_accessor(atags::scalar(), in_(a.derived())),
				n, p, s);
		return p;
	}

	template<class A, typename T>
	inline typename std::enable_if<supports_linear_macc<A>::value,
	index_t>::type
	find_imin(const IEWiseMatrix<A, T>& a)
	{
		const index_t n = a.nelems();
		if (n == 0)
			throw invalid_argument("find_imin: argument a was empty.");

		index_t p;
		T s;
		internal::_find_min(make_vec_accessor(atags::scalar(), in_(a.derived())),
				n, p, s);
		return p;
	}

	template<class A, typename T>
	inline typename std::enable_if<supports_linear_macc<A>::value,
	std::pair<index_t, T> >::type
	find_max(const IEWiseMatrix<A, T>& a)
	{
		const index_t n = a.nelems();
		if (n == 0)
			throw invalid_argument("find_max: argument a was empty.");

		index_t p;
		T s;
		internal::_find_max(make_vec_accessor(atags::scalar(), in_(a.derived())),
				n, p, s);
		return std::make_pair(p, s);
	}

	template<class A, typename T>
	inline typename std::enable_if<supports_linear_macc<A>::value,
	std::pair<index_t, T> >::type
	find_min(const IEWiseMatrix<A, T>& a)
	{
		const index_t n = a.nelems();
		if (n == 0)
			throw invalid_argument("find_min: argument a was empty.");

		index_t p;
		T s;
		internal::_find_min(make_vec_accessor(atags::scalar(), in_(a.derived())),
				n, p, s);
		return std::make_pair(p, s);
	}


	template<typename T, class A, typename TI, class D>
	inline typename std::enable_if<meta::supports_linear_index<D>::value,
	void>::type
	colwise_find_imax(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TI>& idx)
	{
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		if (m == 0)
			throw invalid_argument("colwise_find_imax: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();

		index_t p;
		T s;

		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_max(rd.col(j), m, p, s);
			idx_[j] = static_cast<TI>(p);
		}
	}

	template<typename T, class A, typename TI, class D>
	inline typename std::enable_if<meta::supports_linear_index<D>::value,
	void>::type
	colwise_find_imin(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TI>& idx)
	{
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		if (m == 0)
			throw invalid_argument("colwise_find_imin: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();

		index_t p;
		T s;

		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_min(rd.col(j), m, p, s);
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
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		if (m == 0)
			throw invalid_argument("colwise_find_imax: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();
		R& r_ = r.derived();

		index_t p;
		T s;

		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_max(rd.col(j), m, p, s);
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
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		if (m == 0)
			throw invalid_argument("colwise_find_imin: argument a was empty.");

		LMAT_CHECK_DIMS( a.ncolumns() == idx.nelems() )
		D& idx_ = idx.derived();
		R& r_ = r.derived();

		index_t p;
		T s;

		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			internal::_find_min(rd.col(j), m, p, s);
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
		inline T _nth_elem(IRegularMatrix<A, T>& a, index_t k)
		{
			T *p = a.ptr_data();
			index_t n = a.nelems();

			std::nth_element(p, p+k, p+n);
			return p[k];
		}

		template<typename T, class A>
		inline T _median(IRegularMatrix<A, T>& a)
		{
			T *p = a.ptr_data();
			index_t n = a.nelems();

			T r;
			index_t k = (n >> 1);
			std::nth_element(p, p+k, p+n);

			if (n & 1) // odd number
			{
				r = p[k];
			}
			else  // even number
			{
				T v2 = p[k];
				T v1 = *(std::max_element(p, p+k));
				r = v1 + (v2 - v1) / T(2);
			}

			return r;
		}

	}


	template<class A, typename T>
	inline T nth_element(const IMatrixXpr<A, T>& a, index_t k)
	{
		index_t n = a.nelems();
		if ( k < 0 || k >= n )
			throw invalid_argument("nth_element: the value of k is out of valid range.");

		dense_matrix<T, meta::nrows<A>::value, meta::ncols<A>::value> tmp(a);
		return internal::_nth_elem(tmp, k);
	}


	template<class A, typename T, class D>
	inline typename std::enable_if<meta::supports_linear_index<D>::value,
	void>::type
	colwise_nth_element(const IMatrixXpr<A, T>& a, index_t k, IRegularMatrix<D, T>& r)
	{
		index_t m = a.nrows();
		index_t n = a.ncolumns();
		if ( k < 0 || k >= m )
			throw invalid_argument("colwise_nth_element: the value of k is out of valid range.");

		dense_matrix<T, meta::nrows<A>::value, meta::ncols<A>::value> tmp(a);
		D& r_ = r.derived();
		for (index_t j = 0; j < n; ++j)
		{
			auto cj = tmp.column(j);
			r_[j] = internal::_nth_elem(cj, k);
		}
	}


	template<class A, typename T>
	inline T median(const IMatrixXpr<A, T>& a)
	{
		index_t n = a.nelems();
		if (n == 0)
			throw invalid_argument("median: the input array a was emtpy.");

		dense_matrix<T, meta::nrows<A>::value, meta::ncols<A>::value> tmp(a);
		return internal::_median(tmp);
	}

	template<class A, typename T, class D>
	inline typename std::enable_if<meta::supports_linear_index<D>::value,
	void>::type
	colwise_median(const IMatrixXpr<A, T>& a, IRegularMatrix<D, T>& r)
	{
		if (is_empty(a))
			throw invalid_argument("median: the input array a was emtpy.");

		const index_t n = a.ncolumns();
		LMAT_CHECK_DIMS( n == r.nelems() )

		dense_matrix<T, meta::nrows<A>::value, meta::ncols<A>::value> tmp(a);
		D& r_ = r.derived();
		for (index_t j = 0; j < n; ++j)
		{
			auto cj = tmp.column(j);
			r_[j] = internal::_median(cj);
		}
	}

}

#endif 
