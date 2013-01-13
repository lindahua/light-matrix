/**
 * @file gamma_distr.h
 *
 * @brief Gamma distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_GAMMA_DISTR_H_
#define LIGHTMAT_GAMMA_DISTR_H_

#include "internal/gamma_distr_internal.h"

namespace lmat { namespace random {

	// classes

	template<typename T, typename Method>
	class std_gamma_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		std_gamma_distr(const T& alpha)
		: m_impl(alpha) { }

		LMAT_ENSURE_INLINE
		T alpha() const
		{
			return m_impl.alpha();
		}

		LMAT_ENSURE_INLINE
		T beta() const
		{
			return T(1);
		}

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return alpha();
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return alpha();
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return m_impl(rs);
		}

	private:
		internal::std_gamma_distr_impl<T, Method> m_impl;
	};


	template<typename T, typename Method>
	class gamma_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		gamma_distr(const T& alpha, const T& beta=T(1))
		: m_impl(alpha), m_beta(beta) { }

		LMAT_ENSURE_INLINE
		T alpha() const
		{
			return m_impl.alpha();
		}

		LMAT_ENSURE_INLINE
		T beta() const
		{
			return m_beta;
		}

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return alpha() * m_beta;
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return alpha() * math::sqr(m_beta);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return m_impl(rs) * m_beta;
		}

	private:
		internal::std_gamma_distr_impl<T, Method> m_impl;
		T m_beta;
	};

} }


#endif
