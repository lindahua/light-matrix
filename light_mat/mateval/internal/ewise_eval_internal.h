/**
 * @file ewise_eval_internal.h
 *
 * @brief Internal implementation of element-wise evaluation
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_EWISE_EVAL_INTERNAL_H_
#define LIGHTMAT_EWISE_EVAL_INTERNAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/mateval/vec_accessors.h>
#include <light_mat/mateval/multicol_accessors.h>

namespace lmat { namespace internal {


	/********************************************
	 *
	 *  linear element-wise evaluation
	 *
	 ********************************************/

	template<int N, class Kernel, typename... Accessors>
	inline void linear_ewise_eval_a(const dimension<N>& dim, atags::scalar,
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

	template<int N, typename T, typename SKind, class Kernel, typename... Accessors>
	inline void linear_ewise_eval_a(const dimension<N>& dim, atags::simd<T, SKind>,
			const Kernel& kernel, const Accessors&... accessors)
	{
		const unsigned int W = math::simd_traits<T, SKind>::pack_width;

		const index_t len = dim.value();
		const index_t maj_len = int_div<W>::maj(len);

		if (maj_len)
		{
			pass(accessors.begin_packs()...);

			for (index_t i = 0; i < maj_len; i += W)
			{
				kernel(accessors.pack(i)...);
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


	template<int N, typename U, class Kernel, typename... Wraps>
	inline void linear_ewise_eval(const dimension<N>& dim, U u,
			const Kernel& kernel, const Wraps&... wraps)
	{
		linear_ewise_eval_a(dim, u, kernel, make_vec_accessor(u, wraps)...);
	}


	/********************************************
	 *
	 *  percol element-wise evaluation
	 *
	 ********************************************/

	template<int M, int N, typename U, class Kernel, typename... Accessors>
	inline void percol_ewise_eval_a(const matrix_shape<M, N>& shape, U u,
			const Kernel& kernel, const Accessors&... mcol_accessors)
	{
		const index_t n = shape.ncolumns();
		dimension<M> col_dim(shape.nrows());

		for (index_t j = 0; j < n; ++j)
		{
			linear_ewise_eval_a(col_dim, u, kernel, mcol_accessors.col(j)...);
		}

		pass(mcol_accessors.finalize()...);
	}

	template<int M, int N, typename U, class Kernel, typename... Wraps>
	inline void percol_ewise_eval(const matrix_shape<M, N>& shape, U u,
			const Kernel& kernel, const Wraps&... wraps)
	{
		percol_ewise_eval_a(shape, u, kernel, make_multicol_accessor(u, wraps)...);
	}

} }

#endif
