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

	// traits

	template<typename T, typename Kind>
	struct distr_supp_simd<std_uniform_real_distr<T>, Kind>
	{
		static const bool value = true;
	};

	template<typename T, typename Kind>
	struct distr_supp_simd<uniform_real_distr<T>, Kind>
	{
		static const bool value = true;
	};


	// SIMD generator

	template<typename T, typename Kind>
	class std_uniform_simd_rng
	{
	public:
		typedef math::simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit std_uniform_simd_rng()
		: m_one(T(1)) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		pack_type operator() (RStream& rs) const
		{
			pack_type pk = internal::randbits_to_c1o2(rs.rand_pack(Kind()));
			return pk - m_one;
		}

	private:
		pack_type m_one;
	};


	template<typename T, typename Kind>
	class uniform_simd_rng
	{
		typedef math::simd_pack<T, Kind> pack_type;

		LMAT_ENSURE_INLINE
		explicit uniform_simd_rng(T base, T span)
		: m_base(base), m_span(span) { }

		template<class RStream>
		LMAT_ENSURE_INLINE
		pack_type operator() (RStream& rs) const
		{
			pack_type pk = internal::randbits_to_c1o2(rs.rand_pack(Kind()));
			return m_base + pk * m_span;
		}

	private:
		pack_type m_base;
		pack_type m_span;
	};



	// classes

	template<typename T>
	class std_uniform_real_distr  // U [0, 1)
	{
	public:
		LMAT_ENSURE_INLINE
		T a() const { return T(0); }

		LMAT_ENSURE_INLINE
		T b() const { return T(1); }

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (const RStream& rs)
		{
			return internal::randbits_to_c1o2(rs.rand_u32()) - T(1);
		}
	};


	template<typename T>
	class uniform_real_distr
	{
	public:
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

		template<class RStream>
		LMAT_ENSURE_INLINE
		T operator() (const RStream& rs)
		{
			return m_base + internal::randbits_to_c1o2(rs.rand_u64()) * m_span;
		}

	private:
		T m_a;
		T m_b;
		T m_base;
		T m_span;
	};


} }

#endif
