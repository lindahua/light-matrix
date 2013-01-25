/**
 * @file matrix_find.h
 *
 * @brief Find elements under specific conditions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_FIND_H_
#define LIGHTMAT_MATRIX_FIND_H_

#include <light_mat/matrix/matrix_classes.h>
#include "internal/matrix_find_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  counting
	 *
	 ********************************************/

	template<class A, typename T>
	inline size_t count(const IEWiseMatrix<A, T>& a)
	{
		return internal::count_impl<supports_linear_access<A>::value>::run(a.derived());
	}

	template<class A, typename T, class D, typename TD>
	inline void colwise_count(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TD>& dmat)
	{
		D& d = dmat.derived();
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		LMAT_CHECK_DIMS( d.nelems() == n )
		auto rd = make_multicol_accessor(scalar_(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			d[j] = static_cast<TD>(internal::count_by_reader(m, rd.col(j)));
		}
	}


	/********************************************
	 *
	 *  finding
	 *
	 ********************************************/

	template<class A, typename T, class Visitor>
	inline void findl_f(const IEWiseMatrix<A, T>& a, Visitor vis)
	{
		auto rd = make_vec_accessor(scalar_(), in_(a.derived()));

		const index_t n = a.nelems();
		for (index_t i = 0; i < n; ++i)
		{
			T v = rd.scalar(i);
			if (v) vis(i, v);
		}
	}

	template<class A, typename T, class Visitor>
	inline void find_f(const IEWiseMatrix<A, T>& a, Visitor vis)
	{
		auto rd = make_multicol_accessor(scalar_(), in_(a.derived()));

		const index_t m = a.nrows();
		const index_t n = a.ncolumns();
		for (index_t j = 0; j < n; ++j)
		{
			auto rj = rd.col(j);
			for (index_t i = 0; i < m; ++i)
			{
				T v = rj.scalar(i);
				if (v) vis(i, j, v);
			}
		}
	}


	template<class A, typename T, class VecI>
	inline void findl_to(const IEWiseMatrix<A, T>& a, VecI& veci)
	{
		findl_f(a, [&veci](const index_t& i, const T& ) {
			veci.push_back(i); } );
	}

	template<class A, typename T, class VecI, class VecV>
	inline void findl_to(const IEWiseMatrix<A, T>& a, VecI& veci, VecV& vecv)
	{
		findl_f(a, [&veci, &vecv](const index_t& i, const T& v) {
			veci.push_back(i);
			vecv.push_back(v); } );
	}

	template<class A, typename T, class VecI, class VecJ>
	inline void find_to(const IEWiseMatrix<A, T>& a, VecI& veci, VecJ& vecj)
	{
		find_f(a, [&veci, &vecj](const index_t& i, const index_t& j, const T& ) {
			veci.push_back(i);
			vecj.push_back(j); } );
	}

	template<class A, typename T, class VecI, class VecJ, class VecV>
	inline void find_to(const IEWiseMatrix<A, T>& a, VecI& veci, VecJ& vecj, VecV& vecv)
	{
		find_f(a, [&veci, &vecj, &vecv](const index_t& i, const index_t& j, const T& v) {
			veci.push_back(i);
			vecj.push_back(j);
			vecv.push_back(v); } );
	}



}

#endif /* MATRIX_FIND_H_ */
