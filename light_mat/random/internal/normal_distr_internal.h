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

#include <light_mat/random/uniform_real_distr.h>
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
		typedef T result_type;

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return math::norminv(rand_real<T>::o0c1(rs));
		}
	};

	template<typename T>
	struct normal_distr_impl<T, icdf_>
	{
		typedef T result_type;

		T m_mu;
		T m_sigma;

		LMAT_ENSURE_INLINE
		explicit normal_distr_impl(const T& mu, const T& sigma)
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
			return m_mu + math::norminv(rand_real<T>::o0c1(rs)) * m_sigma;
		}
	};


} } }

#endif
