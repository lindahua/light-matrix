/**
 * @file discrete_distr.h
 *
 * @brief Discrete distribution PRNG
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DISCRETE_DISTR_H_
#define LIGHTMAT_DISCRETE_DISTR_H_

#include <light_mat/random/uniform_int_distr.h>
#include <random>

namespace lmat { namespace random {

	namespace internal
	{
		template<typename TI, typename Method>
		struct discrete_distr_impl;
	}

	template<typename TI, typename Method>
	class discrete_distr  // Uniform over [0, n)
	{
		typedef internal::discrete_distr_impl<TI, Method> impl_t;

	public:
		typedef TI result_type;

		template<typename InputIter>
		LMAT_ENSURE_INLINE
		explicit discrete_distr(InputIter first, InputIter last)
		: m_impl(first, last) { }



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
