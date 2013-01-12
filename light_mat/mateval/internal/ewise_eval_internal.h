/**
 * @file ewise_eval_internal.h
 *
 * @brief Internal implementation of element-wise evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EWISE_EVAL_INTERNAL_H_
#define LIGHTMAT_EWISE_EVAL_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/mateval/vec_accessors.h>
#include <light_mat/mateval/multicol_accessors.h>
#include <light_mat/math/functor_base.h>

namespace lmat { namespace internal {


	/********************************************
	 *
	 *  linear element-wise evaluation
	 *
	 ********************************************/

	template<int N, class Kernel, typename... Accessors>
	inline void linear_ewise_eval(const dimension<N>& dim, atags::scalar,
			const Kernel& kernel, const Accessors&... accessors)
	{
		const index_t len = dim.value();
		for (index_t i = 0; i < len; ++i)
		{
			kernel(accessors.scalar(i)...);
			pass(accessors.done_scalar(i)...);
		}

		pass(accessors.finalize()...);
	}

	template<int N, typename SKind, class Kernel, typename... Accessors>
	inline void linear_ewise_eval(const dimension<N>& dim, atags::simd<SKind>,
			const Kernel& kernel, const Accessors&... accessors)
	{
		static_assert(is_simdizable<Kernel, SKind>::value, "kernel must be simdizable.");

		typedef typename Kernel::value_type T;
		const unsigned int W = simd_traits<T, SKind>::pack_width;

		const index_t len = dim.value();
		const index_t maj_len = static_cast<index_t>(int_div<W>::maj(static_cast<size_t>(len)));

		if (maj_len)
		{
			auto pk_kernel = lmat::simdize_map<Kernel, SKind>::get(kernel);

			pass(accessors.begin_packs()...);

			for (index_t i = 0; i < maj_len; i += W)
			{
				pk_kernel(accessors.pack(i)...);
				pass(accessors.done_pack(i)...);
			}

			pass(accessors.end_packs()...);
		}

		for (index_t i = maj_len; i < len; ++i)
		{
			kernel(accessors.scalar(i)...);
			pass(accessors.done_scalar(i)...);
		}

		pass(accessors.finalize()...);
	}


	/********************************************
	 *
	 *  per-column evaluation
	 *
	 ********************************************/

	template<int N, class VecKernel, typename... MultiColAccessors>
	inline void percol_eval_a(const dimension<N>& col_dim, const index_t n,
			const VecKernel& veckernel, const MultiColAccessors&... accs)
	{
		for (index_t j = 0; j < n; ++j)
		{
			veckernel.apply(col_dim, accs.col(j)...);
		}
	}

} }

#endif
