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

		LMAT_ENSURE_INLINE
		void update(const T& x)
		{
			min_value = math::min(min_value, x);
			max_value = math::max(max_value, x);
		}

		LMAT_ENSURE_INLINE
		void update(const minmax_stat& a)
		{
			min_value = math::min(min_value, a.min_value);
			max_value = math::max(max_value, a.max_value);
		}

		LMAT_ENSURE_INLINE
		T distance() const
		{
			return max_value - min_value;
		}
	};

	template<typename T, typename Kind>
	struct minmax_stat_pack
	{
		typedef simd_pack<T, Kind> pack_t;
		static const unsigned int pack_width = pack_t::pack_width;
		pack_t min_pack;
		pack_t max_pack;

		LMAT_ENSURE_INLINE
		void update(const pack_t& x)
		{
			min_pack = math::min(min_pack, x);
			max_pack = math::max(max_pack, x);
		}

		LMAT_ENSURE_INLINE
		void update(const minmax_stat_pack& a)
		{
			min_pack = math::min(min_pack, a.min_pack);
			max_pack = math::max(max_pack, a.max_pack);
		}
	};


	template<typename T>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> minmax_empty_value()
	{
		T v0 = std::numeric_limits<T>::infinity();
		T v1 = -std::numeric_limits<T>::infinity();
		return minmax_stat<T>({v0, v1});
	}

	template<typename T>
	struct minmax_kernel
	{
		typedef T value_type;
		typedef minmax_stat<T> accumulated_type;

		LMAT_ENSURE_INLINE
		accumulated_type init(const T& x) const
		{
			return accumulated_type({x, x});
		}

		LMAT_ENSURE_INLINE
		void operator()(accumulated_type& a, const T& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		void operator()(accumulated_type& a, const accumulated_type& x) const
		{
			a.update(x);
		}
	};

	template<typename T, typename Kind>
	struct minmax_kernel<simd_pack<T, Kind> >
	{
		typedef simd_pack<T, Kind> value_type;
		typedef minmax_stat_pack<T, Kind> accumulated_type;

		LMAT_ENSURE_INLINE
		accumulated_type init(const value_type& x) const
		{
			return accumulated_type({x, x});
		}

		LMAT_ENSURE_INLINE
		void operator() (accumulated_type& a, const value_type& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		void operator() (accumulated_type& a, const accumulated_type& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		minmax_stat<T> reduce(const accumulated_type& a) const
		{
			T v0 = minimum(a.min_pack);
			T v1 = maximum(a.max_pack);
			return minmax_stat<T>({v0, v1});
		}
	};

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
