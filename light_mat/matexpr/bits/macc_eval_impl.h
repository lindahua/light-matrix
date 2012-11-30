/**
 * @file macc_eval_impl.h
 *
 * @brief Internal implementation of macc-based evaluation
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MACC_EVAL_IMPL_H_
#define LIGHTMAT_MACC_EVAL_IMPL_H_

#include <light_mat/matexpr/matrix_access_base.h>
#include "macc_eval_core.h"

namespace lmat { namespace internal {

	struct macc_eval_linear_scalar
	{
		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static void evaluate(index_t nelems, const SExpr& sexpr, DMat& dmat)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(ct_supports_linear_index<DMat>::value, "DMat must support linear indexing.");
#endif
			typedef typename
					macc_accessor_map<SExpr, linear_macc, scalar_kernel_t>::type
					accessor_t;

			accessor_t accessor(sexpr);
			LinVecRW<DMat> out(dmat);

			typedef macc_vec_copy<scalar_kernel_t, common_ctsize<SExpr, DMat>::value> vec_impl_t;
			vec_impl_t::eval(dmat.nelems(), accessor, out);
		}
	};


	struct macc_eval_percol_scalar
	{
		template<class SExpr, class DMat>
		static void evaluate(index_t m, index_t n, const SExpr& sexpr, DMat& dmat)
		{
			typedef typename
					macc_accessor_map<SExpr, percol_macc, scalar_kernel_t>::type
					accessor_t;

			typedef typename percol_macc_state_map<accessor_t>::type col_state_t;
			typedef typename matrix_traits<DMat>::value_type T;

			typedef percol_to_linear_accessor<accessor_t, scalar_kernel_t, T> acc_wrapper;
			typedef macc_vec_copy<scalar_kernel_t, common_ctrows<SExpr, DMat>::value> vec_impl_t;


			accessor_t accessor(sexpr);

			if (n == 1)
			{
				col_state_t s = accessor.col_state(0);
				T *pd = dmat.ptr_data();
				index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					ContVecRW<T> out(pd);
					vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
				}
				else
				{
					StepVecRW<T> out(pd, rs);
					vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
				}
			}
			else if (m == 1)
			{
				T *pd = dmat.ptr_data();
				index_t cs = dmat.col_stride();

				if (cs == 1)
				{
					for (index_t j = 0; j < n; ++j)
						pd[j] = accessor.get_scalar(0, accessor.col_state(j));
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						pd[j * cs] = accessor.get_scalar(0, accessor.col_state(j));
				}
			}
			else
			{
				T *pd = dmat.ptr_data();
				index_t rs = dmat.row_stride();
				index_t cs = dmat.col_stride();

				if (rs == 1)
				{
					for (index_t j = 0; j < n; ++j, pd += cs)
					{
						col_state_t s = accessor.col_state(j);
						ContVecRW<T> out(pd);
						vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j, pd += cs)
					{
						col_state_t s = accessor.col_state(j);
						StepVecRW<T> out(pd, rs);
						vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
					}
				}
			}
		}
	};

} }

#endif /* MACC_EVAL_IMPL_H_ */
