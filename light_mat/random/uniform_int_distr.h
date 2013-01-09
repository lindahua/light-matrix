/**
 * @file uniform_int_distr.h
 *
 * @brief Uniform integer distribution
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_UNIFORM_INT_DISTR_H_
#define LIGHTMAT_UNIFORM_INT_DISTR_H_

#include <light_mat/random/distr_fwd.h>

namespace lmat { namespace random {

	namespace internal
	{
		template<class RStream, typename T, unsigned int S>
		struct rand_int_helper;

		template<class RStream, typename T>
		struct rand_int_helper<RStream, T, 4>
		{
			LMAT_ENSURE_INLINE
			static T get(RStream& rs)
			{
				return static_cast<T>(rs.rand_u32());
			}
		};

		template<class RStream>
		struct rand_int_helper<RStream, uint32_t, 4>
		{
			LMAT_ENSURE_INLINE
			static uint32_t get(RStream& rs)
			{
				return rs.rand_u32();
			}
		};

		template<class RStream, typename T>
		struct rand_int_helper<RStream, T, 8>
		{
			LMAT_ENSURE_INLINE
			static T get(RStream& rs)
			{
				return static_cast<T>(rs.rand_u64());
			}
		};

		template<class RStream>
		struct rand_int_helper<RStream, uint64_t, 4>
		{
			LMAT_ENSURE_INLINE
			static uint64_t get(RStream& rs)
			{
				return rs.rand_u64();
			}
		};

		template<class RStream, typename T>
		struct rand_int_helper<RStream, T, 1>
		{
			LMAT_ENSURE_INLINE
			static T get(RStream& rs)
			{
				return static_cast<T>(rs.rand_u32() & 0x000000ffU);
			}
		};

		template<class RStream, typename T>
		struct rand_int_helper<RStream, T, 2>
		{
			LMAT_ENSURE_INLINE
			static T get(RStream& rs)
			{
				return static_cast<T>(rs.rand_u32() & 0x0000ffffU);
			}
		};

		template<class RStream, typename T>
		LMAT_ENSURE_INLINE
		inline T get_rand_int(RStream& rs)
		{
			return rand_int_helper<RStream, T, (unsigned int)sizeof(T)>::get(rs);
		}
	}


	template<typename TI>
	class std_uniform_int_distr  // Uniform over [0, b]
	{
	public:
		typedef TI result_type;

	public:
		explicit std_uniform_int_distr(const TI& b)
		: m_span(b+1) { }

		LMAT_ENSURE_INLINE
		TI a() const { return 0; }

		LMAT_ENSURE_INLINE
		TI b() const { return m_span - 1; }

		LMAT_ENSURE_INLINE
		TI span() const { return m_span; }

		LMAT_ENSURE_INLINE
		double p(TI x) const
		{
			return is_nonneg_int(x) && x < m_span ? 1.0 / double(m_span) : 0.0;
		}

		LMAT_ENSURE_INLINE
		double mean() const
		{
			return double(b()) * 0.5;
		}

		LMAT_ENSURE_INLINE
		double var() const
		{
			return (math::sqr(double(m_span)) - 1.0) * (1.0/12);
		}

		template<class RStream>
		LMAT_ENSURE_INLINE
		TI operator() (RStream& rs) const
		{
			return internal::get_rand_int<RStream, TI>(rs) % m_span;
		}

	private:
		TI m_span;
	};


	template<typename TI>
	class uniform_int_distr  // Uniform over [a, b]
	{
	public:
		typedef TI result_type;

	public:
		uniform_int_distr(const TI& a, const TI& b)
		: m_a(a), m_b(b), m_span(b-a+1) { }

		LMAT_ENSURE_INLINE
		TI a() const { return m_a; }

		LMAT_ENSURE_INLINE
		TI b() const { return m_b; }

		LMAT_ENSURE_INLINE
		double p(TI x) const
		{
			return x >= m_a && x <= m_b ? 1.0 / double(m_span) : 0.0;
		}

		LMAT_ENSURE_INLINE
		double mean() const
		{
			return (double(m_a) + double(m_b)) * 0.5;
		}

		LMAT_ENSURE_INLINE
		double var() const
		{
			return (math::sqr(double(m_span)) - 1.0) * (1.0/12);
		}

		LMAT_ENSURE_INLINE
		TI span() const { return m_span; }

		template<class RStream>
		LMAT_ENSURE_INLINE
		TI operator() (RStream& rs) const
		{
			return (internal::get_rand_int<RStream, TI>(rs) % m_span) + m_a;
		}

	private:
		TI m_a;
		TI m_b;
		TI m_span;
	};


} }

#endif
