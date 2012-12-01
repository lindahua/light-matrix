/**
 * @file partial_reduce.h
 *
 * Partial reduction expression and evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PARTIAL_REDUCE_H_
#define LIGHTMAT_PARTIAL_REDUCE_H_

#include <light_mat/matexpr/matexpr_fwd.h>

#include <light_mat/math/reduction_functors.h>
#include <light_mat/matexpr/matrix_arith.h>

#include "bits/partial_reduce_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename Tag, typename Arg_HP, class Arg>
	struct matrix_traits< unary_partial_reduce_expr<Tag, colwise, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = 1;
		static const int ct_num_cols = ct_cols<Arg>::value;
		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef typename reduction_result<Tag, arg_value_type>::type value_type;

		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Tag, typename Arg_HP, class Arg>
	struct matrix_traits< unary_partial_reduce_expr<Tag, rowwise, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value;
		static const int ct_num_cols = 1;
		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef typename reduction_result<Tag, arg_value_type>::type value_type;

		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct matrix_traits< binary_partial_reduce_expr<Tag, colwise, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = 1;
		static const int ct_num_cols = common_ctcols<Arg1, Arg2>::value;
		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename reduction_result<Tag, arg1_value_type, arg2_value_type>::type value_type;

		typedef typename binary_domain<Arg1, Arg2>::type domain;
	};

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct matrix_traits< binary_partial_reduce_expr<Tag, rowwise, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = common_ctrows<Arg1, Arg2>::value;
		static const int ct_num_cols = 1;
		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;
		typedef typename reduction_result<Tag, arg1_value_type, arg2_value_type>::type value_type;

		typedef typename binary_domain<Arg1, Arg2>::type domain;
	};


	/********************************************
	 *
	 *  expression classes
	 *
	 ********************************************/

	namespace internal
	{
		template<class Arg>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_nrows(colwise, const Arg& )
		{
			return 1;
		}

		template<class Arg1, class Arg2>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_nrows(colwise, const Arg1& , const Arg2& )
		{
			return 1;
		}

		template<class Arg>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_ncolumns(colwise, const Arg& a)
		{
			return a.ncolumns();
		}

		template<class Arg1, class Arg2>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_ncolumns(colwise, const Arg1& a1, const Arg2& a2)
		{
			return get_common_ncolumns(a1, a2);
		}

		template<class Arg>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_nrows(rowwise, const Arg& a)
		{
			return a.nrows();
		}

		template<class Arg1, class Arg2>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_nrows(rowwise, const Arg1& a1, const Arg2& a2)
		{
			return get_common_nrows(a1, a2);
		}

		template<class Arg>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_ncolumns(rowwise, const Arg&)
		{
			return 1;
		}

		template<class Arg1, class Arg2>
		LMAT_ENSURE_INLINE
		inline index_t partial_reduc_ncolumns(rowwise, const Arg1&, const Arg2& )
		{
			return 1;
		}
	}


	template<typename Tag, typename AlongDim, typename Arg_HP, class Arg>
	class unary_partial_reduce_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
		unary_partial_reduce_expr<Tag, AlongDim, Arg_HP, Arg>,
		typename matrix_traits< unary_partial_reduce_expr<Tag, AlongDim, Arg_HP, Arg> >::value_type >
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		typedef typename matrix_traits< unary_partial_reduce_expr >::shape_type shape_type;

		LMAT_ENSURE_INLINE
		unary_partial_reduce_expr(const Tag& t,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_tag(t)
		, m_shape( internal::partial_reduc_nrows(AlongDim(), arg_fwd.arg),
				   internal::partial_reduc_ncolumns(AlongDim(), arg_fwd.arg) )
		{ }

		LMAT_ENSURE_INLINE const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

	private:
		Tag m_tag;
		shape_type m_shape;
	};


	template<typename Tag, typename AlongDim, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_partial_reduce_expr
	: public binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IMatrixXpr<
		binary_partial_reduce_expr<Tag, AlongDim, Arg1_HP, Arg1, Arg2_HP, Arg2>,
		typename matrix_traits< binary_partial_reduce_expr<Tag, AlongDim, Arg1_HP, Arg1, Arg2_HP, Arg2> >::value_type >
	{
	public:
		typedef binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;
		typedef typename matrix_traits< binary_partial_reduce_expr >::shape_type shape_type;

		LMAT_ENSURE_INLINE
		binary_partial_reduce_expr(const Tag& t,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		: base_t(arg1_fwd, arg2_fwd), m_tag(t)
		, m_shape( internal::partial_reduc_nrows(AlongDim(), arg1_fwd.arg, arg2_fwd.arg),
				   internal::partial_reduc_ncolumns(AlongDim(), arg1_fwd.arg, arg2_fwd.arg) )
		{ }

		LMAT_ENSURE_INLINE const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

	private:
		Tag m_tag;
		shape_type m_shape;
	};


	/********************************************
	 *
	 *  expression construction
	 *
	 ********************************************/

	// spec classes

	template<typename Tag, typename AlongDim>
	struct par_reduc_t
	{
		Tag tag;

		LMAT_ENSURE_INLINE
		par_reduc_t(const Tag& t) : tag(t) { }
	};

	template<typename Tag, typename AlongDim>
	LMAT_ENSURE_INLINE
	inline par_reduc_t<Tag, AlongDim> par_reduc(const Tag& tag, AlongDim)
	{
		return tag;
	}

	// verifiers

	template<typename Tag, typename AlongDim, class Arg>
	struct unary_expr_verifier<par_reduc_t<Tag, AlongDim>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<typename Tag, typename AlongDim, class Arg1, class Arg2>
	struct binary_expr_verifier<par_reduc_t<Tag, AlongDim>, Arg1, Arg2>
	{
		static const bool value =
				is_mat_xpr<Arg1>::value && is_mat_xpr<Arg2>::value &&
				has_compatible_ct_size<Arg1, Arg2>::value;
	};

	// maps

	template<typename Tag, typename AlongDim, typename Arg_HP, class Arg>
	struct unary_expr_map<par_reduc_t<Tag, AlongDim>, Arg_HP, Arg>
	{
		typedef unary_partial_reduce_expr<Tag, AlongDim, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const par_reduc_t<Tag, AlongDim>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.tag, arg_fwd);
		}
	};

	template<typename Tag, typename AlongDim, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct binary_expr_map<par_reduc_t<Tag, AlongDim>, Arg1_HP, Arg1, Arg2_HP, Arg2>
	{
		typedef binary_partial_reduce_expr<Tag, AlongDim, Arg1_HP, Arg1, Arg2_HP, Arg2> type;

		LMAT_ENSURE_INLINE
		static type get(const par_reduc_t<Tag, AlongDim>& spec,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		{
			return type(spec.tag, arg1_fwd, arg2_fwd);
		}
	};

	// expression construction

	template<typename Tag, typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<par_reduc_t<Tag, colwise>,
		ref_arg_t, Arg
	>::type
	reduce(const Tag& tag, colwise, const IMatrixXpr<Arg, T>& arg)
	{
		return make_expr( par_reduc(tag, colwise()), ref_arg(arg.derived()) );
	}

	template<typename Tag, typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<par_reduc_t<Tag, rowwise>,
		ref_arg_t, Arg
	>::type
	reduce(const Tag& tag, rowwise, const IMatrixXpr<Arg, T>& arg)
	{
		return make_expr( par_reduc(tag, rowwise()), ref_arg(arg.derived()) );
	}

	template<typename Tag, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<par_reduc_t<Tag, colwise>,
		ref_arg_t, Arg1,
		ref_arg_t, Arg2
	>::type
	reduce(  const Tag& tag, colwise,
			const IMatrixXpr<Arg1, T1>& arg1,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr(par_reduc(tag, colwise()),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}

	template<typename Tag, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<par_reduc_t<Tag, rowwise>,
		ref_arg_t, Arg1,
		ref_arg_t, Arg2
	>::type
	reduce(  const Tag& tag, rowwise,
			const IMatrixXpr<Arg1, T1>& arg1,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr(par_reduc(tag, rowwise()),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}


	/********************************************
	 *
	 *  specific expressions
	 *
	 ********************************************/

#define LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( Fname, Tag ) \
		template<typename T, class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<par_reduc_t<Tag, colwise>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, T>& arg, colwise) \
		{ return reduce( Tag(), colwise(), arg ); } \
		template<typename T, class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<par_reduc_t<Tag, rowwise>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, T>& arg, rowwise) \
		{ return reduce( Tag(), rowwise(), arg ); }

#define LMAT_DEFINE_BINARY_PARTIAL_REDUCE_FUNCTION( Fname, Tag ) \
		template<typename T, class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<par_reduc_t<Tag, colwise>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const IMatrixXpr<Arg2, T>& arg2, colwise ) \
		{ return reduce( Tag(), colwise(), arg1, arg2 ); } \
		template<typename T, class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<par_reduc_t<Tag, rowwise>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const IMatrixXpr<Arg2, T>& arg2, rowwise ) \
		{ return reduce( Tag(), rowwise(), arg1, arg2 ); }


	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( sum, sum_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( maximum, maximum_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( minimum, minimum_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( mean, mean_t )

	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( L1norm, L1norm_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( sqL2norm, sqL2norm_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( L2norm, L2norm_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( Linfnorm, Linfnorm_t )

	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( logsum, logsum_t )
	LMAT_DEFINE_UNARY_PARTIAL_REDUCE_FUNCTION( entropy, entropy_t )

	LMAT_DEFINE_BINARY_PARTIAL_REDUCE_FUNCTION( dot, dot_t )
	LMAT_DEFINE_BINARY_PARTIAL_REDUCE_FUNCTION( nrmdot, nrmdot_t )


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<typename Tag, typename AlongDim>
	struct par_reduc_scheme
	{
		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(
				const unary_partial_reduce_expr<Tag, AlongDim, Arg_HP, Arg>& sexpr,
				DMat& dmat)
		{
			internal::_partial_reduce(sexpr.tag(), AlongDim(),
					sexpr.arg(), dmat, scalar_ker());
		}

		template<typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(
				const binary_partial_reduce_expr<Tag, AlongDim, Arg1_HP, Arg1, Arg2_HP, Arg2>& sexpr,
				DMat& dmat)
		{
			internal::_partial_reduce(sexpr.tag(), AlongDim(),
					sexpr.first_arg(), sexpr.second_arg(), dmat, scalar_ker());
		}
	};


	template<typename Tag, typename AlongDim,
		typename Arg_HP, class Arg, class DMat, typename TD>
	LMAT_ENSURE_INLINE
	inline par_reduc_scheme<Tag, AlongDim>
	get_default_eval_scheme(
			const unary_partial_reduce_expr<Tag, AlongDim, Arg_HP, Arg>& sexpr,
			IDenseMatrix<DMat, TD>& dmat)
	{
		return par_reduc_scheme<Tag, AlongDim>();
	}

	template<typename Tag, typename AlongDim,
		typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2, class DMat, typename TD>
	LMAT_ENSURE_INLINE
	inline par_reduc_scheme<Tag, AlongDim>
	get_default_eval_scheme(
			const binary_partial_reduce_expr<Tag, AlongDim, Arg1_HP, Arg1, Arg2_HP, Arg2>& sexpr,
			IDenseMatrix<DMat, TD>& dmat)
	{
		return par_reduc_scheme<Tag, AlongDim>();
	}


}

#endif
