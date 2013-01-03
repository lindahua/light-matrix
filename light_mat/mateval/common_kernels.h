/**
 * @file common_kernels.h
 *
 * @brief commonly used kernels
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_COMMON_KERNELS_H_
#define LIGHTMAT_COMMON_KERNELS_H_

#include <light_mat/mateval/mateval_fwd.h>

#include <light_mat/math/sse_ops.h>
#ifdef LMAT_HAS_AVX
#include <light_mat/math/avx_ops.h>
#endif

namespace lmat
{
	template<typename T>
	struct copy_kernel
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		void operator() (const T& s, T& d) const
		{
			d = s;
		}

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void operator() (const math::simd_pack<T, Kind>& s, math::simd_pack<T, Kind>& d) const
		{
			d = s;
		}
	};

	template<typename Fun>
	struct map_kernel
	{
		typedef typename Fun::value_type value_type;
		Fun fun;

		LMAT_ENSURE_INLINE
		map_kernel(const Fun& f) : fun(f) { }

		template<typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (value_type& out, const Tin&... in) const
		{
			out = fun(in...);
		}

		template<typename Kind, typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<value_type, Kind>& out, const math::simd_pack<Tin, Kind>&... in) const
		{
			out = fun(in...);
		}
	};


	template<typename T>
	struct accum_kernel
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& x) const
		{
			a += x;
		}

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<T, Kind>& a, const math::simd_pack<T, Kind>& x) const
		{
			a += x;
		}
	};

	template<typename T>
	struct accumx_kernel
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& c, const T& x) const
		{
			a += x * c;
		}

		template<typename Kind>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<T, Kind>& a,
				const math::simd_pack<T, Kind>& c, const math::simd_pack<T, Kind>& x) const
		{
			a += x * c;
		}
	};

	template<typename Fun>
	struct accumf_kernel
	{
		typedef typename Fun::value_type value_type;
		Fun fun;

		LMAT_ENSURE_INLINE
		accumf_kernel(const Fun& f) : fun(f) { }

		template<typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (value_type& a, const Tin&... x) const
		{
			a += fun(x...);
		}

		template<typename Kind, typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<value_type, Kind>& a, const math::simd_pack<Tin, Kind>&... in) const
		{
			a += fun(in...);
		}
	};

	template<typename Fun>
	struct accumfx_kernel
	{
		typedef typename Fun::value_type value_type;
		Fun fun;

		LMAT_ENSURE_INLINE
		accumfx_kernel(const Fun& f) : fun(f) { }

		template<typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (value_type& a, const value_type& c, const Tin&... x) const
		{
			a += c * fun(x...);
		}

		template<typename Kind, typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (math::simd_pack<value_type, Kind>& a,
				const math::simd_pack<value_type, Kind>& c, const math::simd_pack<Tin, Kind>&... in) const
		{
			a += c * fun(in...);
		}
	};

}


#endif /* COMMON_KERNELS_H_ */
