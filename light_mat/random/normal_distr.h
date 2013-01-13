/**
 * @file normal_distr.h
 *
 * @brief Normal distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_NORMAL_DISTR_H_
#define LIGHTMAT_NORMAL_DISTR_H_

#include "internal/normal_distr_internal.h"

namespace lmat { namespace random {

	// classes

	template<typename T, typename Method>
	class std_normal_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return T(0);
		}

		LMAT_ENSURE_INLINE
		T stddev() const
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
			return m_impl(rs);
		}

	private:
		internal::std_normal_distr_impl<T, Method> m_impl;
	};


	template<typename T, typename Method>
	class normal_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		explicit normal_distr(T mu=T(0), T sigma=T(1))
		: m_impl(mu, sigma)
		{ }

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return m_impl.mean();
		}

		LMAT_ENSURE_INLINE
		T stddev() const
		{
			return m_impl.stddev();
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return math::sqr(stddev());
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return m_impl(rs);
		}

	private:
		internal::normal_distr_impl<T, Method> m_impl;
	};

} }


namespace lmat
{
	template<typename T, typename Kind>
	struct is_simdizable<random::std_normal_distr<T, random::icdf_>, Kind>
	: public meta::has_simd_support<ftags::norminv_, T, Kind> { };

	template<typename T, typename Kind>
	struct is_simdizable<random::normal_distr<T, random::icdf_>, Kind>
	: public meta::has_simd_support<ftags::norminv_, T, Kind> { };


	template<typename T, typename Kind>
	struct simdize_map< random::std_normal_distr<T, random::icdf_>, Kind >
	{
		typedef random::internal::std_normal_distr_impl<simd_pack<T, Kind>, random::icdf_> type;

		LMAT_ENSURE_INLINE
		static type get(const random::std_normal_distr<T, random::icdf_>& s)
		{
			return type();
		}
	};

	template<typename T, typename Kind>
	struct simdize_map< random::normal_distr<T, random::icdf_>, Kind >
	{
		typedef random::internal::normal_distr_impl<simd_pack<T, Kind>, random::icdf_> type;

		LMAT_ENSURE_INLINE
		static type get(const random::normal_distr<T, random::icdf_>& s)
		{
			return type(s.mean(), s.stddev());
		}
	};

}



#endif
