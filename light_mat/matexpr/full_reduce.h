/**
 * @file full_reduce.h
 *
 * Full matrix reduction expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REDUCE_H_
#define LIGHTMAT_MATRIX_REDUCE_H_

#include "bits/full_reduce_internal.h"

#define LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION( Tag, FunName ) \
	template<typename T, class Xpr> \
	LMAT_ENSURE_INLINE \
	inline typename reductor_traits<Tag, LMAT_TYPELIST_1(T) >::result_type \
	FunName(const IMatrixXpr<Xpr, T>& a) { return reduce(Tag(), a); }

#define LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION( Tag, FunName ) \
	template<typename T, class Xpr1, class Xpr2> \
	LMAT_ENSURE_INLINE \
	inline typename reductor_traits<Tag, LMAT_TYPELIST_2(T, T) >::result_type \
	FunName(const IMatrixXpr<Xpr1, T>& a1, const IMatrixXpr<Xpr2, T>& a2) { return reduce(Tag(), a1, a2); }

#define LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( Name ) \
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION( Name##_t, Name )

#define LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION_( Name ) \
	LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION( Name##_t, Name )


namespace lmat
{
	/********************************************
	 *
	 *  generic reduction functions
	 *
	 ********************************************/

	template<class Tag, typename T, class Arg, class Acc>
	LMAT_ENSURE_INLINE
	typename reductor_traits<Tag, LMAT_TYPELIST_1(T) >::result_type
	reduce(Tag, const IMatrixXpr<Arg, T>& a, macc_policy<Acc, scalar_ker>)
	{
		typedef reductor_traits<Tag, LMAT_TYPELIST_1(T) > traits;
		typedef typename traits::result_type RT;

		RT r;

		if (!is_empty(a))
		{
			typedef typename op_fun<
					typename traits::combine_fun_tag,
					scalar_ker,
					LMAT_TYPELIST_2(RT, RT) >::type combfun_t;

			typedef typename traits::post_fun_tag post_tag_t;

			typedef internal::full_reduc_helper<Tag, LMAT_TYPELIST_1(Arg) > helper_t;
			typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;

			typename helper_t::term_expr_t term_expr =
					helper_t::get_term_expr( tfwd_t(a.derived()) );

			typedef macc_policy<Acc, scalar_ker> policy_t;


			internal::full_reduce_(Tag(), term_expr, combfun_t(), post_tag_t(), r, policy_t());
		}
		else
		{
			r = traits::empty();
		}

		return r;
	}

	template<class Tag, typename T1, class Arg1, typename T2, class Arg2, class Acc>
	LMAT_ENSURE_INLINE
	typename reductor_traits<Tag, LMAT_TYPELIST_2(T1, T2) >::result_type
	reduce(Tag, const IMatrixXpr<Arg1, T1>& a1, const IMatrixXpr<Arg2, T2>& a2, macc_policy<Acc, scalar_ker>)
	{
		typedef reductor_traits<Tag, LMAT_TYPELIST_2(T1, T2) > traits;
		typedef typename traits::result_type RT;

		RT r;

		const index_t m = common_nrows(a1, a2);
		const index_t n = common_ncols(a1, a2);

		if (m > 0 && n > 0)
		{
			typedef typename op_fun<
					typename traits::combine_fun_tag,
					scalar_ker,
					LMAT_TYPELIST_2(RT, RT) >::type combfun_t;

			typedef typename traits::post_fun_tag post_tag_t;

			typedef internal::full_reduc_helper<Tag, LMAT_TYPELIST_2(Arg1, Arg2) > helper_t;
			typedef tied_forwarder< LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2>) > tfwd_t;

			typename helper_t::term_expr_t term_expr =
					helper_t::get_term_expr( tfwd_t(a1.derived(), a2.derived()) );

			typedef macc_policy<Acc, scalar_ker> policy_t;


			internal::full_reduce_(Tag(), term_expr, combfun_t(), post_tag_t(), r, policy_t());
		}
		else
		{
			r = traits::empty();
		}

		return r;
	}



	template<class Tag, typename T, class Arg>
	LMAT_ENSURE_INLINE
	typename reductor_traits<Tag, LMAT_TYPELIST_1(T) >::result_type
	reduce(Tag, const IMatrixXpr<Arg, T>& a)
	{
		typedef reductor_traits<Tag, LMAT_TYPELIST_1(T) > traits;
		typedef typename traits::result_type RT;

		RT r;

		if (!is_empty(a))
		{
			typedef typename op_fun<
					typename traits::combine_fun_tag,
					scalar_ker,
					LMAT_TYPELIST_2(RT, RT) >::type combfun_t;

			typedef typename traits::post_fun_tag post_tag_t;

			typedef internal::full_reduc_helper<Tag, LMAT_TYPELIST_1(Arg) > helper_t;
			typedef typename helper_t::term_expr_t term_expr_t;
			typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;

			term_expr_t texpr = helper_t::get_term_expr( tfwd_t(a.derived()) );

			if (helper_t::decide_use_linear(texpr))
			{
				typedef macc_policy<linear_macc, scalar_ker> policy_t;
				internal::full_reduce_(Tag(), texpr, combfun_t(), post_tag_t(), r, policy_t());
			}
			else
			{
				typedef macc_policy<percol_macc, scalar_ker> policy_t;
				internal::full_reduce_(Tag(), texpr, combfun_t(), post_tag_t(), r, policy_t());
			}
		}
		else
		{
			r = traits::empty();
		}

		return r;
	}


	template<class Tag, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	typename reductor_traits<Tag, LMAT_TYPELIST_2(T1, T2) >::result_type
	reduce(Tag, const IMatrixXpr<Arg1, T1>& a1, const IMatrixXpr<Arg2, T2>& a2)
	{
		typedef reductor_traits<Tag, LMAT_TYPELIST_2(T1, T2) > traits;
		typedef typename traits::result_type RT;

		RT r;

		const index_t m = common_nrows(a1, a2);
		const index_t n = common_ncols(a1, a2);

		if (m > 0 && n > 0)
		{
			typedef typename op_fun<
					typename traits::combine_fun_tag,
					scalar_ker,
					LMAT_TYPELIST_2(RT, RT) >::type combfun_t;

			typedef typename traits::post_fun_tag post_tag_t;

			typedef internal::full_reduc_helper<Tag, LMAT_TYPELIST_2(Arg1, Arg2) > helper_t;
			typedef typename helper_t::term_expr_t term_expr_t;
			typedef tied_forwarder< LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2>) > tfwd_t;

			term_expr_t texpr = helper_t::get_term_expr( tfwd_t(a1.derived(), a2.derived()) );

			if (helper_t::decide_use_linear(texpr))
			{
				typedef macc_policy<linear_macc, scalar_ker> policy_t;
				internal::full_reduce_(Tag(), texpr, combfun_t(), post_tag_t(), r, policy_t());
			}
			else
			{
				typedef macc_policy<percol_macc, scalar_ker> policy_t;
				internal::full_reduce_(Tag(), texpr, combfun_t(), post_tag_t(), r, policy_t());
			}
		}
		else
		{
			r = traits::empty();
		}

		return r;
	}



	/********************************************
	 *
	 *  specific reduction functions
	 *
	 ********************************************/

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( sum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( maximum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( minimum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( mean )

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( L1norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( sqL2norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( L2norm )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( Linfnorm )

	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( logsum )
	LMAT_DEFINE_UNARY_FULL_REDUCE_FUNCTION_( entropy )

	LMAT_DEFINE_BINARY_FULL_REDUCE_FUNCTION_( dot )
}

#endif /* MATRIX_REDUC_EXPR_H_ */


