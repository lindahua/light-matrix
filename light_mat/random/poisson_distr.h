/**
 * @file poisson_distr.h
 *
 * @brief Poisson distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_POISSON_DISTR_H_
#define LIGHTMAT_POISSON_DISTR_H_

#include <light_mat/random/exponential_distr.h>

namespace lmat { namespace random {


	/********************************************
	 *
	 *  internal implementation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename TI, typename Method> struct poisson_distr_impl;

		template<typename TI>
		struct poisson_distr_impl<TI, naive_>
		{
		public:
			poisson_distr_impl(double mu)
			: m_mu(mu) { }

			LMAT_ENSURE_INLINE
			double mean() const
			{
				return m_mu;
			}

			template<class RStream>
			LMAT_ENSURE_INLINE
			TI operator() (RStream& rs) const
			{
				TI x = 0;

				double s = 0.0;
				do
				{
					s += m_egen(rs);
					++ x;

				} while ( s < m_mu );

				return x - 1;
			}

		private:
			double m_mu;
			std_exponential_distr<double> m_egen;
		};
	}


	/********************************************
	 *
	 *  distribution class
	 *
	 ********************************************/

	template<typename TI, typename Method>
	class poisson_distr
	{
		typedef internal::poisson_distr_impl<TI, Method> impl_t;

	public:
		typedef TI result_type;

		LMAT_ENSURE_INLINE
		explicit poisson_distr(double mu = 1.0)
		: m_impl(mu)
		{ }

		double p(TI x) const
		{
			double v;

			if (is_nonneg_int(x))
			{
				const double lam = mean();
				v = math::exp(-lam);

				if (x > 0)
				{
					for (TI i = 1; i <= x; ++i)
						v *= (lam / double(i));
				}
			}
			else
			{
				v = 0.0;
			}

			return v;
		}

		LMAT_ENSURE_INLINE
		double mean() const
		{
			return m_impl.mean();
		}

		LMAT_ENSURE_INLINE
		double var() const
		{
			return m_impl.mean();
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
