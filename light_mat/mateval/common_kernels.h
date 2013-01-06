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
#include <light_mat/math/simd_ops.h>

/************************************************
 *
 *  Macros to define common kernels
 *
 ************************************************/

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
	};

	LMAT_DEF_SIMD_SUPPORT( copy_kernel )

	template<typename Fun>
	struct map_kernel
	{
		typedef typename Fun::result_type value_type;
		Fun fun;

		LMAT_ENSURE_INLINE
		map_kernel(const Fun& f) : fun(f) { }

		template<typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (value_type& out, const Tin&... in) const
		{
			out = fun(in...);
		}
	};

	template<typename Fun>
	struct is_simdizable<map_kernel<Fun> >
	{
		static const bool value = is_simdizable<Fun>::value;
	};

	template<typename Fun, typename Kind>
	struct simdize_map<map_kernel<Fun>, Kind>
	{
		typedef map_kernel<typename simdize_map<Fun, Kind>::type> type;

		LMAT_ENSURE_INLINE
		static type get(const map_kernel<Fun>& kernel)
		{
			return type(simdize_map<Fun, Kind>::get(kernel.fun));
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
	};

	LMAT_DEF_SIMD_SUPPORT( accum_kernel )

	template<typename T>
	struct accumx_kernel
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& c, const T& x) const
		{
			a += x * c;
		}
	};

	LMAT_DEF_SIMD_SUPPORT( accumx_kernel )

}


#endif /* COMMON_KERNELS_H_ */
