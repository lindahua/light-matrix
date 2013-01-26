/**
 * @file svec_internal.h
 *
 * Internal implementation of svec
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SVEC_INTERNAL_H_
#define LIGHTMAT_SVEC_INTERNAL_H_

#include <light_mat/simd/simd.h>
#include "simd_sarith.h"

namespace lmat { namespace internal {

	template<typename T, typename Kind, unsigned int PW, index_t Len>
	class svec_impl;


	/********************************************
	 *
	 *  pack-width = 2, vec-length = 2
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 2, 2>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x);
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1)
		{
			m_pk.set(x0, x1);
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_u(p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_a(p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_u(p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_a(p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


	/********************************************
	 *
	 *  pack-width = 2, vec-length = 3
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 2, 3>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE
		svec_impl(const pack_t& pk0, const pack_t& pk1)
		: m_pk0(pk0), m_pk1(pk1) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk0.reset();
			m_pk1.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk0.set(x);
			m_pk1.set(x, T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2)
		{
			m_pk0.set(x0, x1);
			m_pk1.set(x2, T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk0.load_u(p);
			m_pk1.load_part(siz_<1>(), p + 2);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk0.load_a(p);
			m_pk1.load_part(siz_<1>(), p + 2);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk0.store_u(p);
			m_pk1.store_part(siz_<1>(), p + 2);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk0.store_a(p);
			m_pk1.store_part(siz_<1>(), p + 2);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk0 + x.m_pk0, sadd(m_pk1, x.m_pk1));
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk0 - x.m_pk0, ssub(m_pk1, x.m_pk1));
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk0 * x.m_pk0, smul(m_pk1, x.m_pk1));
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk0 += x.m_pk0;
			m_pk1 = sadd(m_pk1, x.m_pk1);
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk0 -= x.m_pk0;
			m_pk1 = ssub(m_pk1, x.m_pk1);
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk0 *= x.m_pk0;
			m_pk1 = smul(m_pk1, x.m_pk1);
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum(sadd(m_pk0, m_pk1));
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum(sadd(m_pk0 * x.m_pk0, m_pk1 * x.m_pk1));
		}

	private:
		pack_t m_pk0;
		pack_t m_pk1;
	};


	/********************************************
	 *
	 *  pack-width = 2, vec-length = 4
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 2, 4>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE
		svec_impl(const pack_t& pk0, const pack_t& pk1)
		: m_pk0(pk0), m_pk1(pk1) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk0.reset();
			m_pk1.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk0.set(x);
			m_pk1.set(x);
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2, const T& x3)
		{
			m_pk0.set(x0, x1);
			m_pk1.set(x2, x3);
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk0.load_u(p);
			m_pk1.load_u(p + 2);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk0.load_a(p);
			m_pk1.load_a(p + 2);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk0.store_u(p);
			m_pk1.store_u(p + 2);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk0.store_a(p);
			m_pk1.store_a(p + 2);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk0 + x.m_pk0, m_pk1 + x.m_pk1);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk0 - x.m_pk0, m_pk1 - x.m_pk1);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk0 * x.m_pk0, m_pk1 * x.m_pk1);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk0 += x.m_pk0;
			m_pk1 += x.m_pk1;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk0 -= x.m_pk0;
			m_pk1 -= x.m_pk1;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk0 *= x.m_pk0;
			m_pk1 *= x.m_pk1;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum(m_pk0 + m_pk1);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum(m_pk0 * x.m_pk0 + m_pk1 * x.m_pk1);
		}

	private:
		pack_t m_pk0;
		pack_t m_pk1;
	};



	/********************************************
	 *
	 *  pack-width = 4, vec-length = 2
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 4, 2>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x, x, T(0), T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1)
		{
			m_pk.set(x0, x1, T(0), T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_part(siz_<2>(), p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum_2(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum_2(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


	/********************************************
	 *
	 *  pack-width = 4, vec-length = 3
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 4, 3>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x, x, x, T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2)
		{
			m_pk.set(x0, x1, x2, T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_part(siz_<3>(), p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum_3(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum_3(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


	/********************************************
	 *
	 *  pack-width = 4, vec-length = 4
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 4, 4>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x);
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2, const T& x3)
		{
			m_pk.set(x0, x1, x2, x3);
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_u(p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_a(p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_u(p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_a(p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};



	/********************************************
	 *
	 *  pack-width = 8, vec-length = 2
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 8, 2>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x, x, T(0), T(0), T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1)
		{
			m_pk.set(x0, x1, T(0), T(0), T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_part(siz_<2>(), p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_part(siz_<2>(), p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum_2(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum_2(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


	/********************************************
	 *
	 *  pack-width = 8, vec-length = 3
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 8, 3>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x, x, x, T(0), T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2)
		{
			m_pk.set(x0, x1, x2, T(0), T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_part(siz_<3>(), p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_part(siz_<3>(), p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum_3(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum_3(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


	/********************************************
	 *
	 *  pack-width = 8, vec-length = 4
	 *
	 ********************************************/

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 8, 4>
	{
		typedef simd_pack<T, Kind> pack_t;

		LMAT_ENSURE_INLINE svec_impl(const pack_t& pk)
		: m_pk(pk) { }

	public:
		LMAT_ENSURE_INLINE svec_impl() { }

		// set, load & store

		LMAT_ENSURE_INLINE void reset()
		{
			m_pk.reset();
		}

		LMAT_ENSURE_INLINE void set(const T& x)
		{
			m_pk.set(x, x, x, x, T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void set(const T& x0, const T& x1, const T& x2, const T& x3)
		{
			m_pk.set(x0, x1, x2, x3, T(0), T(0), T(0), T(0));
		}

		LMAT_ENSURE_INLINE void load_u(const T *p)
		{
			m_pk.load_part(siz_<4>(), p);
		}

		LMAT_ENSURE_INLINE void load_a(const T *p)
		{
			m_pk.load_part(siz_<4>(), p);
		}

		LMAT_ENSURE_INLINE void store_u(T *p) const
		{
			m_pk.store_part(siz_<4>(), p);
		}

		LMAT_ENSURE_INLINE void store_a(T *p) const
		{
			m_pk.store_part(siz_<4>(), p);
		}

		// arithmetics

		LMAT_ENSURE_INLINE svec_impl operator + (const svec_impl& x) const
		{
			return svec_impl(m_pk + x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator - (const svec_impl& x) const
		{
			return svec_impl(m_pk - x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl operator * (const svec_impl& x) const
		{
			return svec_impl(m_pk * x.m_pk);
		}

		LMAT_ENSURE_INLINE svec_impl& operator += (const svec_impl& x)
		{
			m_pk += x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator -= (const svec_impl& x)
		{
			m_pk -= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE svec_impl& operator *= (const svec_impl& x)
		{
			m_pk *= x.m_pk;
			return *this;
		}

		LMAT_ENSURE_INLINE T _sum() const
		{
			return sum_4(m_pk);
		}

		LMAT_ENSURE_INLINE T _dot(const svec_impl& x) const
		{
			return sum_4(m_pk * x.m_pk);
		}

	private:
		pack_t m_pk;
	};


} }

#endif 
