/**
 * @file full_reduce_internal.h
 *
 * @brief Internal implementation of full matrix reduction
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_FULL_REDUCE_INTERNAL_H_
#define LIGHTMAT_FULL_REDUCE_INTERNAL_H_

#include <light_mat/math/reductors.h>
#include <light_mat/matexpr/dense_accessors.h>
#include <light_mat/matexpr/matrix_ewise_eval.h>

namespace lmat { namespace internal {

	template<class Tag, class ArgList>
	struct full_reduc_helper
	{
		static const int M = meta::common_nrows<ArgList>::value;
		static const int N = meta::common_ncols<ArgList>::value;

		typedef typename meta::transform_<
				meta::value_type_of, ArgList>::type arg_value_types;

		typedef reductor_traits<Tag, arg_value_types> traits;
		typedef typename traits::result_type result_type;

		typedef typename meta::transform_<
				meta::to_cref_arg, ArgList>::type qargs;

		typedef typename traits::term_fun_tag term_fun_tag;
		typedef ewise_expr<term_fun_tag, qargs> term_expr_t;

		LMAT_ENSURE_INLINE
		static term_expr_t get_term_expr( const tied_forwarder<qargs>& tfwd )
		{
			return term_expr_t(term_fun_tag(), tfwd);
		}

		LMAT_ENSURE_INLINE
		static bool decide_use_linear( const term_expr_t& sexpr )
		{
			const index_t m = sexpr.nrows();

			int macc_linear_cost = macc_cost<term_expr_t, linear_macc, scalar_ker>::value;
			int macc_percol_cost = macc_cost<term_expr_t, percol_macc, scalar_ker>::value;

			if (m <= MACC_SHORT_PERCOL_COST)
				macc_percol_cost += MACC_SHORT_PERCOL_COST;

			return macc_linear_cost <= macc_percol_cost;
		}
	};


	template<typename RT, int M, int N>
	struct full_reduce_impl
	{
		static const int CTSize = M * N;

		template<class CombFun, class Acc>
		LMAT_ENSURE_INLINE
		static RT fold(index_t m, index_t n, const CombFun& combine, const Acc& acc,
				macc_policy<linear_macc, scalar_ker>)
		{
			if (CTSize > 0)
			{
				RT r = acc.get_scalar(0);
				for (index_t i = 1; i < CTSize; ++i)
					r = combine(r, acc.get_scalar(i));
				return r;
			}
			else
			{
				const index_t len = m * n;

				RT r = acc.get_scalar(0);
				for (index_t i = 1; i < len; ++i)
					r = combine(r, acc.get_scalar(i));
				return r;
			}
		}


		template<class CombFun, class Acc>
		static RT fold(index_t m, index_t n, const CombFun& combine, const Acc& acc,
				macc_policy<percol_macc, scalar_ker>)
		{
			typedef typename percol_macc_state_map<Acc>::type col_state_t;

			RT r;

			if (M == 1 || m == 1)
			{
				col_state_t s0(acc.col_state(0));
				r = acc.get_scalar(0, s0);

				for (index_t j = 1; j < n; ++j)
				{
					col_state_t sj(acc.col_state(j));
					r = combine(r, acc.get_scalar(0, sj));
				}
			}
			else
			{
				// first column

				col_state_t s0(acc.col_state(0));
				r = acc.get_scalar(0, s0);

				if (M > 0)
				{
					for (index_t i = 1; i < M; ++i)
						r = combine(r, acc.get_scalar(i, s0));

					for (index_t j = 1; j < n; ++j)
					{
						col_state_t sj(acc.col_state(j));

						RT rj = acc.get_scalar(0, sj);
						for (index_t i = 1; i < M; ++i)
							rj = combine(rj, acc.get_scalar(i, sj));

						r = combine(r, rj);
					}
				}
				else
				{
					for (index_t i = 1; i < m; ++i)
						r = combine(r, acc.get_scalar(i, s0));

					for (index_t j = 1; j < n; ++j)
					{
						col_state_t sj(acc.col_state(j));

						RT rj = acc.get_scalar(0, sj);
						for (index_t i = 1; i < m; ++i)
							rj = combine(rj, acc.get_scalar(i, sj));

						r = combine(r, rj);
					}
				}

			}

			return r;
		}


		template<typename PostTag>
		LMAT_ENSURE_INLINE
		static RT post(PostTag, const RT& folded_val, scalar_ker, index_t, index_t)
		{
			typename op_fun<PostTag, scalar_ker, LMAT_TYPELIST_1(RT) >::type pf;
			return pf(folded_val);
		}

		LMAT_ENSURE_INLINE
		static RT post(id_t, const RT& folded_val, scalar_ker, index_t, index_t)
		{
			return folded_val;
		}

		LMAT_ENSURE_INLINE
		static RT post(divide_by_dim, const RT& folded_val, scalar_ker, index_t m, index_t n)
		{
			index_t len = m * n;
			if (len == 1)
			{
				return folded_val;
			}
			else
			{
				return folded_val / static_cast<RT>(m * n);
			}
		}

	};


	template<class Tag, class TermExpr, class CombFun, class PostFunTag, typename RT, class Policy>
	void full_reduce_(Tag,
			const TermExpr& texpr, const CombFun& combine, const PostFunTag& post_tag,
			RT& r, const Policy& policy)
	{

		typedef full_reduce_impl<RT,
				meta::nrows<TermExpr>::value, meta::ncols<TermExpr>::value> impl_t;

		const index_t m = texpr.nrows();
		const index_t n = texpr.ncolumns();

		typename macc_accessor_map<TermExpr, Policy>::type acc(texpr);

		typedef typename Policy::kernel_kind kernel_t;
		r = impl_t::post(post_tag, impl_t::fold(m, n, combine, acc, policy),
				kernel_t(), m, n);
	}


} }

#endif /* FULL_REDUCE_INTERNAL_H_ */
