/**
 * @file uniform_real_distr.h
 *
 * @brief Real-valued uniform distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_UNIFORM_REAL_DISTR_H_
#define LIGHTMAT_UNIFORM_REAL_DISTR_H_

#include "internal/uniform_real_internal.h"
#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif

namespace lmat { namespace random {

	// classes

	template<typename T>
	class std_uniform_real_distr  // U [0, 1)
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		T a() const { return T(0); }

		LMAT_ENSURE_INLINE
		T b() const { return T(1); }

		LMAT_ENSURE_INLINE
		T span() const { return T(1); }

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return T(0.5);
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return T(1.0/12);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return internal::randreal_helper<T>::c1o2(rs) - T(1);
		}
	};


	template<typename T>
	class uniform_real_distr
	{
	public:
		typedef T result_type;

		LMAT_ENSURE_INLINE
		explicit uniform_real_distr(T a_, T b_)
		: m_a(a_), m_b(b_)
		, m_base(T(2) * a_ - b_)
		, m_span(b_ - a_)
		{ }

		LMAT_ENSURE_INLINE
		T a() const { return m_a; }

		LMAT_ENSURE_INLINE
		T b() const { return m_b; }

		LMAT_ENSURE_INLINE
		T span() const { return m_span; }

		LMAT_ENSURE_INLINE
		T mean() const
		{
			return (m_a + m_b) * T(0.5);
		}

		LMAT_ENSURE_INLINE
		T var() const
		{
			return math::sqr(m_span) * T(1.0/12);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (RStream& rs) const
		{
			return m_base + internal::randreal_helper<T>::c1o2(rs) * m_span;
		}

		LMAT_ENSURE_INLINE
		T _base() const { return m_base; }

	private:
		T m_a;
		T m_b;
		T m_base;
		T m_span;
	};


	// SIMD generator

	template<typename T, typename Kind>
	class std_uniform_real_simd
	{
	public:
		typedef math::simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		explicit std_uniform_real_simd()
		: m_one(T(1)) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return pk - m_one;
		}

	private:
		result_type m_one;
	};


	template<typename T, typename Kind>
	class uniform_real_simd
	{
	public:
		typedef math::simd_pack<T, Kind> result_type;

		LMAT_ENSURE_INLINE
		explicit uniform_real_simd(T base, T span)
		: m_base(base), m_span(span) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		result_type operator() (RStream& rs) const
		{
			result_type pk = internal::randreal_helper<T>::c1o2(rs, Kind());
			return m_base + pk * m_span;
		}

	private:
		result_type m_base;
		result_type m_span;
	};

} }


namespace lmat
{

	LMAT_DECL_SIMDIZABLE_ON_REAL( random::std_uniform_real_distr )
	LMAT_DECL_SIMDIZABLE_ON_REAL( random::uniform_real_distr )

	template<typename T, typename Kind>
	struct simdize_map< random::std_uniform_real_distr<T>, Kind >
	{
		typedef random::std_uniform_real_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::std_uniform_real_distr<T>& s)
		{
			return type();
		}
	};

	template<typename T, typename Kind>
	struct simdize_map< random::uniform_real_distr<T>, Kind >
	{
		typedef random::uniform_real_simd<T, Kind> type;

		LMAT_ENSURE_INLINE
		static type get(const random::uniform_real_distr<T>& s)
		{
			return type(s._base(), s.span());
		}
	};

}


#endif
