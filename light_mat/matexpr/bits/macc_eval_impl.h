/**
 * @file macc_eval_impl.h
 *
 * @brief Internal implementation of macc-based evaluation
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MACC_EVAL_IMPL_H_
#define LIGHTMAT_MACC_EVAL_IMPL_H_

#include <light_mat/matrix/matrix_concepts.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{
	struct linear_macc;
	struct percol_macc;

	template<class SExpr, class AccCate, class KerCate>
	struct macc_accessor_map;

	template<class Accessor>
	struct percol_macc_state_map;

	namespace internal
	{

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

				for (index_t i = 0; i < nelems; ++i)
				{
					dmat[i] = accessor.get_scalar(i);
				}
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

				accessor_t accessor(sexpr);

				if (n == 1)
				{
					col_state_t s = accessor.col_state(0);
					T *pd = dmat.ptr_data();
					index_t rs = dmat.row_stride();

					if (rs == 1)
					{
						for (index_t i = 0; i < m; ++i)
							pd[i] = accessor.get_scalar(i, s);
					}
					else
					{
						for (index_t i = 0; i < m; ++i)
							pd[i * rs] = accessor.get_scalar(i, s);
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
							for (index_t i = 0; i < m; ++i)
								pd[i] = accessor.get_scalar(i, s);
						}
					}
					else
					{
						for (index_t j = 0; j < n; ++j, pd += cs)
						{
							col_state_t s = accessor.col_state(j);
							for (index_t i = 0; i < m; ++i)
								pd[i * rs] = accessor.get_scalar(i, s);
						}
					}
				}
			}
		};

	}

}

#endif /* MACC_EVAL_IMPL_H_ */
