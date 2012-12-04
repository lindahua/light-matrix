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

	template<typename Acc, typename Ker>
	struct macc_eval_impl;

	template<typename Ker>
	struct macc_eval_impl<linear_macc, Ker>
	{
		template<class SExpr, class DMat>
		static void evaluate(index_t nelems, const SExpr& sexpr, DMat& dmat)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(meta::supports_linear_index<DMat>::value, "DMat must support linear indexing.");
#endif
			typedef typename
					macc_accessor_map<SExpr, linear_macc, scalar_ker>::type
					accessor_t;

			typedef typename meta::common_value_type<SExpr, DMat>::type T;

			typedef LMAT_TYPELIST_2( SExpr, DMat ) Lst;

			accessor_t accessor(sexpr);
			T *pd = dmat.ptr_data();

			if (meta::is_continuous<DMat>::value)
			{
				ContVecRW<Ker, T> out(pd);
				typedef macc_vec_copy<Ker, meta::common_nelems<Lst>::value> impl_t;
				impl_t::eval(dmat.nelems(), accessor, out);
			}
			else if (dmat.ncolumns() == 1)
			{
				typedef macc_vec_copy<Ker, meta::common_nrows<Lst>::value> impl_t;
				const index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					ContVecRW<Ker, T> out(pd);
					impl_t::eval(dmat.nrows(), accessor, out);
				}
				else
				{
					StepVecRW<Ker, T> out(pd, rs);
					impl_t::eval(dmat.nrows(), accessor, out);
				}
			}
			else // nrows == 1
			{
				typedef macc_vec_copy<Ker, meta::common_ncols<Lst>::value> impl_t;
				const index_t cs = dmat.col_stride();

				if (cs == 1)
				{
					ContVecRW<Ker, T> out(pd);
					impl_t::eval(dmat.ncolumns(), accessor, out);
				}
				else
				{
					StepVecRW<Ker, T> out(pd, cs);
					impl_t::eval(dmat.ncolumns(), accessor, out);
				}
			}
		}
	};


	template<typename Ker>
	struct macc_eval_impl<percol_macc, Ker>
	{
		template<class SExpr, class DMat>
		static void evaluate(index_t m, index_t n, const SExpr& sexpr, DMat& dmat)
		{
			typedef typename
					macc_accessor_map<SExpr, percol_macc, scalar_ker>::type
					accessor_t;

			typedef typename percol_macc_state_map<accessor_t>::type col_state_t;
			typedef typename matrix_traits<DMat>::value_type T;

			typedef LMAT_TYPELIST_2( SExpr, DMat ) Lst;

			typedef percol_to_linear_accessor<accessor_t, scalar_ker, T> acc_wrapper;
			typedef macc_vec_copy<scalar_ker, meta::common_nrows<Lst>::value> vec_impl_t;

			accessor_t accessor(sexpr);

			if (n == 1)
			{
				col_state_t s = accessor.col_state(0);
				T *pd = dmat.ptr_data();
				index_t rs = dmat.row_stride();

				if (rs == 1)
				{
					ContVecRW<Ker, T> out(pd);
					vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
				}
				else
				{
					StepVecRW<Ker, T> out(pd, rs);
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
						ContVecRW<Ker, T> out(pd);
						vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
					}
				}
				else
				{
					for (index_t j = 0; j < n; ++j, pd += cs)
					{
						col_state_t s = accessor.col_state(j);
						StepVecRW<Ker, T> out(pd, rs);
						vec_impl_t::eval(m, acc_wrapper(accessor, s), out);
					}
				}
			}
		}
	};

} }

#endif /* MACC_EVAL_IMPL_H_ */
