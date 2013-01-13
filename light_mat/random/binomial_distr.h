/**
 * @file binomial_distr.h
 *
 * Binomial distribution
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BINOMIAL_DISTR_H_
#define LIGHTMAT_BINOMIAL_DISTR_H_

#include <light_mat/random/bernoulli_distr.h>


namespace lmat { namespace random {


	/********************************************
	 *
	 *  internal implementation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename TI>
		inline double _nchoosek(TI n, TI k)
		{
			if (k > n-k) k = n-k;

			double x = 1.0;
			double y = 1.0;

			TI b = n-k;

			for (TI i = 1; i <= k; ++i)
			{
				x *= double(b+i);
				y *= double(i);
			}

			return x / y;
		}

		template<typename TI, typename Method> struct binomial_distr_impl;

		template<typename TI>
		struct binomial_distr_impl<TI, naive_>
		{
		public:
			binomial_distr_impl(TI t, double p)
			: m_t(t), m_bernoulli(p) { }

			LMAT_ENSURE_INLINE
			TI t() const
			{
				return m_t;
			}

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
				for (TI i = 0; i < m_t; ++i)
				{
					if (m_bernoulli(rs)) ++x;
				}
				return x;
			}

		private:
			TI m_t;
			bernoulli_distr m_bernoulli;
		};
	}

	/********************************************
	 *
	 *  distribution class
	 *
	 ********************************************/

	template<typename TI, typename Method>
	class binomial_distr
	{
		typedef internal::binomial_distr_impl<TI, Method> impl_t;

	public:
		typedef TI result_type;

		LMAT_ENSURE_INLINE
		explicit binomial_distr(const TI& t, const double& p = 0.5)
		: m_impl(t, p)
		{ }

		LMAT_ENSURE_INLINE
		TI t() const
		{
			return m_impl.t();
		}

		LMAT_ENSURE_INLINE
		double p() const
		{
			return m_impl.p();
		}

		double p(TI k) const
		{
			if ( is_nonneg_int(k) && k <= t() )
			{
				double c = internal::_nchoosek(t(), k);
				TI r = t() - k;

				double x = 1.0;
				double y = 1.0;

				double pp = p();
				double pn = 1.0 - p();

				for (TI i = 0; i < r; ++i) x *= pn;
				for (TI i = 0; i < k; ++i) y *= pp;

				return c * x * y;
			}
			else return 0.0;
		}

		LMAT_ENSURE_INLINE
		double mean() const
		{
			return p() * double(t());
		}

		LMAT_ENSURE_INLINE
		double var() const
		{
			double _p = p();
			return _p * (1.0 - _p) * double(t());
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
