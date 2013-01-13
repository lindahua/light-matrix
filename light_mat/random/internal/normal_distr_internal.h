/**
 * @file normal_distr_internal.h
 *
 * @brief Internal implementation of normal distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_NORMAL_DISTR_INTERNAL_H_
#define LIGHTMAT_NORMAL_DISTR_INTERNAL_H_

#include "uniform_real_internal.h"
#include <light_mat/math/math_special.h>
#include <light_mat/math/simd_math.h>


namespace lmat { namespace random { namespace internal {

	// forward declarations

	template<typename T, typename Method>
	struct std_normal_distr_impl;

	template<typename T, typename Method>
	struct normal_distr_impl;

	template<typename T, typename Kind, typename Method>
	struct std_normal_distr_simd_impl;

	template<typename T, typename Kind, typename Method>
	struct normal_distr_simd_impl;


	/********************************************
	 *
	 *  ICDF implementation
	 *
	 ********************************************/

	template<typename T>
	struct std_normal_distr_impl<T, icdf_>
	{
		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			T u = internal::randreal_helper<T>::c1o2(rs) - T(1);
			return math::norminv(u);
		}
	};

	template<typename T>
	struct normal_distr_impl<T, icdf_>
	{
		T m_mu;
		T m_sigma;

		LMAT_ENSURE_INLINE
		explicit normal_distr_impl(T mu, T sigma)
		: m_mu(mu), m_sigma(sigma)
		{ }

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return m_mu;
		}

		LMAT_ENSURE_INLINE
		T stddev() const
		{
			return m_sigma;
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			T u = internal::randreal_helper<T>::c1o2(rs) - T(1);
			return m_mu + math::norminv(u) * m_sigma;
		}
	};


	template<typename T, typename Kind>
	struct std_normal_distr_simd_impl<T, Kind, icdf_>
	{
		typedef simd_pack<T, Kind> result_type;

		result_type m_one;

		LMAT_ENSURE_INLINE
		explicit std_normal_distr_simd_impl()
		: m_one(T(1)) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return math::norminv(pk - m_one);
		}
	};


	template<typename T, typename Kind>
	struct normal_distr_simd_impl<T, Kind, icdf_>
	{
		typedef simd_pack<T, Kind> result_type;

		result_type m_one;
		result_type m_mu;
		result_type m_sigma;

		LMAT_ENSURE_INLINE
		explicit normal_distr_simd_impl(T mu, T sigma)
		: m_one(T(1)), m_mu(mu), m_sigma(sigma) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return m_mu + math::norminv(pk - m_one) * m_sigma;
		}

	};


} } }

#endif
