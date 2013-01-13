/**
 * @file exponential_distr.h
 *
 * @brief Exponential distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EXPONENTIAL_DISTR_H_
#define LIGHTMAT_EXPONENTIAL_DISTR_H_

#include <light_mat/random/uniform_real_distr.h>
#include <light_mat/math/simd_math.h>

namespace lmat { namespace random {

	// classes

	template<typename T>
	class std_exponential_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T lambda() const
		{
			return T(1);
		}

		LMAT_ENSURE_INLINE
		T beta() const
		{
			return T(1);
		}

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return T(1);
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return T(1);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return - math::log(rand_real<T>::o0c1(rs));
		}
	};


	template<typename T>
	class exponential_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		explicit exponential_distr(T lambda=1.0)
		: m_lambda(lambda), m_beta(math::rcp(lambda))
		{ }

		LMAT_ENSURE_INLINE
		T lambda() const
		{
			return m_lambda;
		}

		LMAT_ENSURE_INLINE
		T beta() const
		{
			return m_beta;
		}

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return m_beta;
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return math::sqr(m_beta);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return (-m_beta) * math::log(rand_real<T>::o0c1(rs));
		}

	private:
		T m_lambda;
		T m_beta;
	};


	// SIMD generators

	template<typename T, typename Kind>
	class std_exponential_distr_simd
	{
	public:
		typedef simd_pack<T, Kind> result_type;

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			return - math::log(rand_real<result_type>::o0c1(rs));
		}
	};


	template<typename T, typename Kind>
	class exponential_distr_simd
	{
	public:
		typedef simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		explicit exponential_distr_simd(const T& beta)
		: m_neg_beta(-beta)
		{ }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			return m_neg_beta * math::log(rand_real<result_type>::o0c1(rs));
		}

	private:
		T m_neg_beta;
	};

} }


namespace lmat
{

	template<typename T, typename Kind>
	struct is_simdizable<random::std_exponential_distr<T>, Kind>
	: public meta::has_simd_support<ftags::log_, T, Kind> { };

	template<typename T, typename Kind>
	struct is_simdizable<random::exponential_distr<T>, Kind>
	: public meta::has_simd_support<ftags::log_, T, Kind> { };

	template<typename T, typename Kind>
	struct simdize_map< random::std_exponential_distr<T>, Kind >
	{
		typedef random::std_exponential_distr_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::std_exponential_distr<T>& s)
		{
			return type();
		}
	};

	template<typename T, typename Kind>
	struct simdize_map< random::exponential_distr<T>, Kind >
	{
		typedef random::exponential_distr_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::exponential_distr<T>& s)
		{
			return type(s.beta());
		}
	};

}


#endif

