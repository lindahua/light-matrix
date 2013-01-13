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

#include "internal/uniform_real_internal.h"
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
			T u = T(2) - internal::randreal_helper<T>::c1o2(rs);
			return - math::log(u);
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
			T u = T(2) - internal::randreal_helper<T>::c1o2(rs);
			return (-m_beta) * math::log(u);
		}

	private:
		T m_lambda;
		T m_beta;
	};


	// SIMD generator

	template<typename T, typename Kind>
	class std_exponential_simd
	{
	public:
		typedef simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		explicit std_exponential_simd()
		: m_two(T(2)) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return - math::log(m_two - pk);
		}

	private:
		result_type m_two;
	};


	template<typename T, typename Kind>
	class exponential_simd
	{
	public:
		typedef simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		explicit exponential_simd(T beta)
		: m_two(T(2)), m_neg_beta(-beta) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return math::log(m_two - pk) * m_neg_beta;
		}

	private:
		result_type m_two;
		result_type m_neg_beta;
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
		typedef random::std_exponential_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::std_exponential_distr<T>& s)
		{
			return type();
		}
	};

	template<typename T, typename Kind>
	struct simdize_map< random::exponential_distr<T>, Kind >
	{
		typedef random::exponential_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::exponential_distr<T>& s)
		{
			return type(s.beta());
		}
	};

}


#endif

