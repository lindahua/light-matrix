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

/************************************************
 *
 *  Macros to define common kernels
 *
 ************************************************/

#define LMAT_pkargt ::lmat::math::simd_pack<T, Kind>&

#define LMAT_kerarg1(P) P##1
#define LMAT_kerarg2(P) P##1, P##2
#define LMAT_kerarg3(P) P##1, P##2, P##3
#define LMAT_kerarg4(P) P##1, P##2, P##3, P##4
#define LMAT_kerarg5(P) P##1, P##2, P##3, P##4, P##5
#define LMAT_kerarg6(P) P##1, P##2, P##3, P##4, P##5, P##6
#define LMAT_kerarg7(P) P##1, P##2, P##3, P##4, P##5, P##6, P##7
#define LMAT_kerarg8(P) P##1, P##2, P##3, P##4, P##5, P##6, P##7, P##8
#define LMAT_kerarg9(P) P##1, P##2, P##3, P##4, P##5, P##6, P##7, P##8, P##9

#define LMAT_simd_pkt lmat::math::simd_pack<T, Kind>

#define LMAT_DEF_GKERNEL(Kernel, Fun, NIn, NOut) \
	template<typename T> \
	struct Kernel { \
		typedef T value_type; \
		LMAT_ENSURE_INLINE \
		void operator() (LMAT_kerarg##NIn(const T& i), LMAT_kerarg##NOut(T& o)) const \
		{ Fun(LMAT_kerarg##NIn(i), LMAT_kerarg##NOut(o)); } \
		template<typename Kind> \
		LMAT_ENSURE_INLINE \
		void operator() (LMAT_kerarg##NIn(const LMAT_simd_pkt& i), LMAT_kerarg##NOut(LMAT_simd_pkt& o)) const \
		{ Fun(LMAT_kerarg##NIn(i), LMAT_kerarg##NOut(o)); } };


namespace lmat
{

	/************************************************
	 *
	 *  specific kernels
	 *
	 ************************************************/

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
