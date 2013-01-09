/**
 * @file geometric_distr.h
 *
 * 
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_GEOMETRIC_DISTR_H_
#define LIGHTMAT_GEOMETRIC_DISTR_H_

#include <light_mat/random/bernoulli_distr.h>

namespace lmat { namespace random {

	/********************************************
	 *
	 *  internal implementation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename TI, typename Method> struct geometric_distr_impl;

		template<typename TI>
		struct geometric_distr_impl<TI, naive_>
		{
		public:
			geometric_distr_impl(TI t, double p)
			: m_bernoulli(p) { }

			LMAT_ENSURE_INLINE
			double p() const
			{
				return m_bernoulli.p();
			}

			template<class RStream>
			LMAT_ENSURE_INLINE
			TI operator() (RStream& rs) const
			{
				TI x(0);
				while (!m_bernoulli(rs)) ++x;
				return x;
			}

		private:
			bernoulli_distr m_bernoulli;
		};
	}

	/********************************************
	 *
	 *  distribution class
	 *
	 ********************************************/

	template<typename TI, typename Method>
	class geometric_distr
	{
		typedef internal::geometric_distr_impl<TI, Method> impl_t;

	public:
		LMAT_ENSURE_INLINE
		explicit geometric_distr(const double& p_ = 0.5)
		: m_impl(p_)
		{ }

		LMAT_ENSURE_INLINE
		double t() const
		{
			return m_impl.t();
		}

		LMAT_ENSURE_INLINE
		double p() const
		{
			return m_impl.p();
		}

		LMAT_ENSURE_INLINE
		double mean() const
		{
			return math::rcp(p()) - 1.0;
		}

		LMAT_ENSURE_INLINE
		double var() const
		{
			double r = math::rcp(p());
			return math::sqr(r) - r;
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		TI operator() (RStream& rs) const
		{
			return m_impl(rs);
		}

	private:
		impl_t m_impl;
	};


} }

#endif 
