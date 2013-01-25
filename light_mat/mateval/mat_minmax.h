/**
 * @file mat_minmax.h
 *
 * @brief min-max reduction of matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_MINMAX_H_
#define LIGHTMAT_MAT_MINMAX_H_

#include "internal/mat_reduce_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  minmax statistics
	 *
	 ********************************************/

	template<typename T>
	struct minmax_stat
	{
		T min_value;
		T max_value;

		minmax_stat() { }

		minmax_stat(const T& x)
		: min_value(x), max_value(x) { }

		void update(const T& x)
		{
			min_value = math::min(min_value, x);
			max_value = math::max(max_value, x);
		}

		void update(const minmax_stat& s)
		{
			min_value = math::min(min_value, s.min_value);
			max_value = math::max(max_value, s.max_value);
		}
	};

	template<typename T>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> minmax_empty_value()
	{
		minmax_stat<T> s;
		s.min_value = std::numeric_limits<T>::infinity();
		s.max_value = - std::numeric_limits<T>::infinity();
		return s;
	}

	template<typename T, typename Kind>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> reduce_impl(const minmax_stat<simd_pack<T, Kind> >& s)
	{
		minmax_stat<T> r;
		r.min_value = minimum(s.min_value);
		r.max_value = maximum(s.max_value);
		return r;
	}

	LMAT_DEFINE_AGGREG_SIMD_FOLDKERNEL(minmax_stat, minmax_kernel, 1)

	LMAT_DEF_SIMD_SUPPORT( minmax_kernel )


	/********************************************
	 *
	 *  minmax reduction functions
	 *
	 ********************************************/

	template<typename T, class A>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> minmax(const IEWiseMatrix<A, T>& a)
	{
		return a.nelems() > 0 ?
				fold(minmax_kernel<T>())(a.shape(), in_(a)) :
				minmax_empty_value<T>();
	}


	template<typename T, class A, class DMat1, class DMat2>
	inline void colwise_minmax(const IEWiseMatrix<A, T>& a,
			IRegularMatrix<DMat1, T>& dmat_min,
			IRegularMatrix<DMat2, T>& dmat_max)
	{
		dimension<meta::nrows<A>::value> col_dim(a.nrows());
		const index_t n = a.ncolumns();
		LMAT_CHECK_DIMS( n == dmat_min.nelems() && n == dmat_max.nelems() )

		auto g = make_colwise_fold_getter(minmax_kernel<T>(), a.shape(), a);

		DMat1& d1 = dmat_min.derived();
		DMat2& d2 = dmat_max.derived();

		for (index_t j = 0; j < n; ++j)
		{
			auto r = g[j];
			d1[j] = r.min_value;
			d2[j] = r.max_value;
		}
	}


}

#endif /* MAT_MINMAX_H_ */
