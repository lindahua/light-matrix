/**
 * @file common_kernels.h
 *
 * @brief commonly used kernels
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_COMMON_KERNELS_H_
#define LIGHTMAT_COMMON_KERNELS_H_

#include <light_mat/mateval/mateval_fwd.h>


namespace lmat
{
	struct copy_kernel
	{
		template<typename T>
		LMAT_ENSURE_INLINE
		void operator() (const T& s, T& d) const
		{
			d = s;
		}
	};

	template<typename Fun>
	struct map_kernel
	{
		Fun fun;

		LMAT_ENSURE_INLINE
		map_kernel(const Fun& f) : fun(f) { }

		template<typename Tout, typename... Tin>
		LMAT_ENSURE_INLINE
		void operator() (Tout& out, const Tin&... in) const
		{
			out = fun(in...);
		}
	};


	struct accum_kernel
	{
		template<typename T>
		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& x) const
		{
			a += x;
		}
	};

	struct accumx_kernel
	{
		template<typename T>
		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& c, const T& x) const
		{
			a += x * c;
		}
	};

	template<typename Fun>
	struct accumf_kernel
	{
		Fun fun;

		LMAT_ENSURE_INLINE
		accumf_kernel(const Fun& f) : fun(f) { }

		template<typename T, typename... A>
		LMAT_ENSURE_INLINE
		void operator() (T& a, const A&... x) const
		{
			a += fun(x...);
		}
	};

	template<typename Fun>
	struct accumfx_kernel
	{
		Fun fun;

		LMAT_ENSURE_INLINE
		accumfx_kernel(const Fun& f) : fun(f) { }

		template<typename T, typename... A>
		LMAT_ENSURE_INLINE
		void operator() (T& a, const T& c, const A&... x) const
		{
			a += c * fun(x...);
		}
	};

}


#endif /* COMMON_KERNELS_H_ */
