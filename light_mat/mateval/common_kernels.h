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

}


#endif /* COMMON_KERNELS_H_ */
