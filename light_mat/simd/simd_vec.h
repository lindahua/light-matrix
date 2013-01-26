/**
 * @file simd_vec.h
 *
 * @brief SIMD short vector classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SIMD_VEC_H_
#define LIGHTMAT_SIMD_VEC_H_

#include "internal/svec_internal.h"

namespace lmat
{
	template<typename T, typename Kind, index_t Len>
	class simd_vec
	: public internal::svec_impl<T, Kind, simd_traits<T, Kind>::pack_width, Len>
	{
	private:
		typedef internal::svec_impl<T, Kind, simd_traits<T, Kind>::pack_width, Len> base_t;

	public:
		typedef T value_type;
		typedef Kind simd_kind;

		LMAT_ENSURE_INLINE
		simd_vec() { }

		LMAT_ENSURE_INLINE
		simd_vec(const base_t& base) : base_t(base) { }

		LMAT_ENSURE_INLINE
		simd_vec(const T& v)
		{
			this->set(v);
		}

		LMAT_ENSURE_INLINE
		simd_vec(const T* src)
		{
			this->load_u(src);
		}

		LMAT_ENSURE_INLINE
		index_t length() const
		{
			return Len;
		}
	};


	template<typename T, typename Kind, index_t Len>
	LMAT_ENSURE_INLINE
	inline T sum(const simd_vec<T, Kind, Len>& vec)
	{
		return vec._sum();
	}

	template<typename T, typename Kind, index_t Len>
	LMAT_ENSURE_INLINE
	inline T dot(const simd_vec<T, Kind, Len>& vx, const simd_vec<T, Kind, Len>& vy)
	{
		return vx._dot(vy);
	}

}

#endif
