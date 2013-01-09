/**
 * @file bernoulli_distr.h
 *
 * @brief Bernoulli distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_BERNOULLI_DISTR_H_
#define LIGHTMAT_BERNOULLI_DISTR_H_

#include <light_mat/random/distr_fwd.h>

namespace lmat { namespace random {

	class std_bernoulli_distr
	{
	public:
		LMAT_ENSURE_INLINE
		std_bernoulli_distr() { }

		LMAT_ENSURE_INLINE
		double p() const { return 0.5; }

		LMAT_ENSURE_INLINE
		double mean() const { return 0.5; }

		LMAT_ENSURE_INLINE
		double var() const { return 0.25; }

		template<class RStream>
		LMAT_ENSURE_INLINE
		bool operator() (RStream& rs) const
		{
			return static_cast<bool>(rs.rand_u32() & 1);
		}
	};


	class bernoulli_distr
	{
	public:
		LMAT_ENSURE_INLINE
		explicit bernoulli_distr(double p=0.5)
		: m_p(p)
		{
			// 4294967296.0 = double(0xffffffff) + 1.0
			// when p = 1.0, it ends up 0xffffffff
			m_thres = static_cast<uint32_t>(4294967296.0 * p);
		}

		LMAT_ENSURE_INLINE
		double p() const { return m_p; }

		LMAT_ENSURE_INLINE
		double mean() const { return m_p; }

		LMAT_ENSURE_INLINE
		double var() const { return m_p * (1.0 - m_p); }

		template<class RStream>
		LMAT_ENSURE_INLINE
		bool operator() (RStream& rs) const
		{
			return rs.rand_u32() < m_thres;
		}

	private:
		double m_p;
		uint32_t m_thres;
	};

} }


#endif
