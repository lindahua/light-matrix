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

namespace lmat { namespace internal {

	template<typename T, typename Kind, unsigned int PW, index_t Len>
	class svec_impl;

	template<typename T, typename Kind>
	class svec_impl<T, Kind, 2, 2>
	{
	public:
		void set_zeros()
		{
			m_pk.reset();
		}

		void set(const T& x)
		{
			m_pk.set(x);
		}

		void load_u(const T *p)
		{
			m_pk.load_u(p);
		}

		void load_a(const T *p)
		{
			m_pk.load_a(p);
		}

		void store_u(T *p) const
		{
			m_pk.store_u(p);
		}

		void store_a(T *p) const
		{
			m_pk.store_a(p);
		}

	private:
		simd_pack<T, Kind> m_pk;
	};

} }

#endif 
