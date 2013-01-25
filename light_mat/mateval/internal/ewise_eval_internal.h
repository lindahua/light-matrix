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

	template<index_t Len, class Kernel, typename... Accessors>
	inline void _linear_ewise_eval(
			const dimension<Len>& dim, scalar_,
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

	template<index_t Len, typename SKind, class Kernel, typename... Accessors>
	inline void _linear_ewise_eval(
			const dimension<Len>& dim, simd_<SKind>,
			const Kernel& kernel, const Accessors&... accessors)
	{
		static_assert(is_simdizable<Kernel, SKind>::value, "kernel must be simdizable.");

		typedef typename Kernel::value_type T;
		const unsigned int W = simd_traits<T, SKind>::pack_width;
		const index_t W_ = static_cast<index_t>(W);

		const index_t len = dim.value();

		if (len >= W_)
		{
			const unsigned int W2 = W * 2;
			const index_t W2_ = static_cast<index_t>(W2);

			const size_t npacks = int_div<W>::quo(static_cast<size_t>(len));
			auto pk_kernel = lmat::simdize_map<Kernel, SKind>::get(kernel);

			index_t maj_len;

			if (npacks > 1)
			{
				maj_len = static_cast<index_t>(int_div<W2>::maj(static_cast<size_t>(len)));
				pass(accessors.begin_packs()...);

				for (index_t i = 0; i < maj_len; i += W2_)
				{
					pk_kernel(accessors.pack(i)...);
					pass(accessors.done_pack(i)...);
					pk_kernel(accessors.pack(i + W_)...);
					pass(accessors.done_pack(i + W_)...);
				}

				if (npacks & 1)
				{
					pk_kernel(accessors.pack(maj_len)...);
					pass(accessors.done_pack(maj_len)...);
					maj_len += W_;
				}

				pass(accessors.end_packs()...);

			}
			else // npacks == 1
			{
				maj_len = W_;
				pass(accessors.begin_packs()...);
				pk_kernel(accessors.pack(0)...);
				pass(accessors.done_pack(0)...);
				pass(accessors.end_packs()...);
			}

			for (index_t i = maj_len; i < len; ++i)
			{
				kernel(accessors.scalar(i)...);
				pass(accessors.done_scalar(i)...);
			}
		}
		else
		{
			for (index_t i = 0; i < len; ++i)
			{
				kernel(accessors.scalar(i)...);
				pass(accessors.done_scalar(i)...);
			}
		}

		pass(accessors.finalize()...);
	}


	/********************************************
	 *
	 *  per-column evaluation
	 *
	 ********************************************/

	template<index_t CM, index_t CN, typename U, class Kernel, typename... MultiColAccessors>
	inline void _percol_ewise_eval(
			const matrix_shape<CM, CN>& shape, U,
			const Kernel& kernel, const MultiColAccessors&... accessors)
	{
		dimension<CM> coldim(shape.nrows());
		const index_t n = shape.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			_linear_ewise_eval(coldim, U(), kernel, accessors.col(j)...);
		}
	}



} }

#endif
