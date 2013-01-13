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
	struct minmax_folder
	{
		typedef minmax_stat<T> value_type;

		LMAT_ENSURE_INLINE
		value_type init(const T& x) const
		{
			return value_type({x, x});
		}

		LMAT_ENSURE_INLINE
		void fold(value_type& a, const T& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		void fold(value_type& a, const value_type& x) const
		{
			a.update(x);
		}
	};


	template<typename T, typename Kind>
	struct minmax_pack_folder
	{
		typedef minmax_stat_pack<T, Kind> value_type;
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE
		value_type init(const simd_pack<T, Kind>& x) const
		{
			return value_type({x, x});
		}

		LMAT_ENSURE_INLINE
		void fold(value_type& a, const pack_t& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		void fold(value_type& a, const value_type& x) const
		{
			a.update(x);
		}

		LMAT_ENSURE_INLINE
		minmax_stat<T> reduce(const value_type& a) const
		{
			T v0 = minimum(a.min_pack);
			T v1 = maximum(a.max_pack);
			return minmax_stat<T>({v0, v1});
		}
	};

	LMAT_DECL_SIMDIZABLE_ON_REAL( minmax_folder )

	template<typename T, typename Kind>
	struct simdize_map<minmax_folder<T>, Kind>
	{
		typedef minmax_pack_folder<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const minmax_folder<T>& )
		{
			return type();
		}
	};


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
				internal::_full_reduce(a.shape(), minmax_folder<T>(), a.derived()) :
				minmax_empty_value<T>();
	}

	template<typename T, class A, typename U>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> minmax(const IEWiseMatrix<A, T>& a, linear_macc<U> policy)
	{
		return a.nelems() > 0 ?
				internal::_full_reduce(a.shape(), minmax_folder<T>(), a.derived(), policy) :
				minmax_empty_value<T>();
	}

	template<typename T, class A, typename U>
	LMAT_ENSURE_INLINE
	inline minmax_stat<T> minmax(const IEWiseMatrix<A, T>& a, percol_macc<U> policy)
	{
		return a.nelems() > 0 ?
				internal::_full_reduce(a.shape(), minmax_folder<T>(), a.derived(), policy) :
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

		typedef minmax_folder<T> folder_t;
		typedef typename internal::colwise_reduc_policy<folder_t, A, DMat1>::atag atag;

		auto fker = fold(folder_t(), atag());
		auto rd = make_multicol_accessor(atag(), in_(a.derived()));

		DMat1& dmat_min_ = dmat_min.derived();
		DMat2& dmat_max_ = dmat_max.derived();

		for (index_t j = 0; j < n; ++j)
		{
			auto s = fker.apply(col_dim, rd.col(j));
			dmat_min_[j] = s.min_value;
			dmat_max_[j] = s.max_value;
		}
	}


}

#endif /* MAT_MINMAX_H_ */
