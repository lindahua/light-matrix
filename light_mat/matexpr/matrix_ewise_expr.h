/*
 * @file matrix_ewise_expr.h
 *
 * Generic element-wise matrix expression
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EWISE_EXPR_H_
#define LIGHTMAT_MATRIX_EWISE_EXPR_H_

#include <light_mat/matexpr/matexpr_fwd.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expression class traits
	 *
	 ********************************************/

	template<typename Op, typename Arg_HP, class Arg>
	struct matrix_traits<unary_ewise_expr<Op, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value;
		static const int ct_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type arg_value_type;

		typedef typename unary_op_result<Op, arg_value_type>::type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Op, typename Arg1_HP, class Arg1, class Arg2_HP, class Arg2>
	struct matrix_traits<binary_ewise_expr<Op, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = common_ctrows<Arg1, Arg2>::value;
		static const int ct_num_cols = common_ctcols<Arg1, Arg2>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;

		typedef typename binary_op_result<Op, arg1_value_type, arg2_value_type>::type value_type;
		typedef typename binary_domain<Arg1, Arg2>::type domain;
	};

	template<typename Op, typename T1, class Arg2_HP, class Arg2>
	struct matrix_traits<binary_fix1st_ewise_expr<Op, T1, Arg2_HP, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg2>::value;
		static const int ct_num_cols = ct_cols<Arg2>::value;

		static const bool is_readonly = true;

		typedef T1 arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;

		typedef typename binary_op_result<Op, arg1_value_type, arg2_value_type>::type value_type;
		typedef typename matrix_traits<Arg2>::domain domain;
	};

	template<typename Op, typename Arg1_HP, class Arg1, typename T2>
	struct matrix_traits<binary_fix2nd_ewise_expr<Op, Arg1_HP, Arg1, T2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg1>::value;
		static const int ct_num_cols = ct_cols<Arg1>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef T2 arg2_value_type;

		typedef typename binary_op_result<Op, arg1_value_type, arg2_value_type>::type value_type;
		typedef typename matrix_traits<Arg1>::domain domain;
	};


	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<typename Op, typename Arg_HP, class Arg>
	class unary_ewise_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<unary_ewise_expr<Op, Arg_HP, Arg>,
		typename matrix_traits<unary_ewise_expr<Op, Arg_HP, Arg> >::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_unary_op<Op>::value, "Op must be a unary operation");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;

		LMAT_ENSURE_INLINE
		unary_ewise_expr(const Op& op,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_op(op) { }

		LMAT_ENSURE_INLINE const Op& op() const
		{
			return m_op;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		Op m_op;
	};


	template<typename Op, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_expr
	: public binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IMatrixXpr<binary_ewise_expr<Op, Arg1_HP, Arg1, Arg2_HP, Arg2>,
		typename matrix_traits<binary_ewise_expr<Op, Arg1_HP, Arg1, Arg2_HP, Arg2> >::value_type >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_op<Op>::value, "Op must be a binary operation.");
		static_assert(is_mat_xpr<Arg1>::value, "Arg1 must be a matrix expression class.");
		static_assert(is_mat_xpr<Arg2>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Op& op,
				const arg_forwarder<Arg1_HP, Arg1>& arg1,
				const arg_forwarder<Arg2_HP, Arg2>& arg2)
		: base_t(arg1, arg2), m_op(op)
		{
			LMAT_CHECK_SAME_SHAPE(this->first_arg(), this->second_arg())
		}

		LMAT_ENSURE_INLINE const Op& op() const
		{
			return m_op;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->first_arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->first_arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->first_arg().ncolumns();
		}

	private:
		Op m_op;
	};


	template<typename Op, typename T1, typename Arg2_HP, class Arg2>
	class binary_fix1st_ewise_expr
	: public unary_expr_base<Arg2_HP, Arg2>
	, public IMatrixXpr<binary_fix1st_ewise_expr<Op, T1, Arg2_HP, Arg2>,
		typename matrix_traits<binary_fix1st_ewise_expr<Op, T1, Arg2_HP, Arg2> >::value_type >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_op<Op>::value, "Op must be a binary operation.");
		static_assert(is_mat_xpr<Arg2>::value, "Arg2 must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg2_HP, Arg2> base_t;

		LMAT_ENSURE_INLINE
		binary_fix1st_ewise_expr(const Op& op,
				const T1& arg1v,
				const arg_forwarder<Arg2_HP, Arg2>& arg2)
		: base_t(arg2), m_op(op), m_arg1v(arg1v)
		{
		}

		LMAT_ENSURE_INLINE const Op& op() const
		{
			return m_op;
		}

		LMAT_ENSURE_INLINE const T1& arg1_value() const
		{
			return m_arg1v;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		Op m_op;
		const T1 m_arg1v;
	};


	template<typename Op, typename Arg1_HP, class Arg1, typename T2>
	class binary_fix2nd_ewise_expr
	: public unary_expr_base<Arg1_HP, Arg1>
	, public IMatrixXpr<binary_fix2nd_ewise_expr<Op, Arg1_HP, Arg1, T2>,
		typename matrix_traits<binary_fix2nd_ewise_expr<Op, Arg1_HP, Arg1, T2> >::value_type >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_op<Op>::value, "Op must be a binary operation.");
		static_assert(is_mat_xpr<Arg1>::value, "Arg1 must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg1_HP, Arg1> base_t;

		LMAT_ENSURE_INLINE
		binary_fix2nd_ewise_expr(const Op& op,
				const arg_forwarder<Arg1_HP, Arg1>& arg1,
				const T2& arg2v)
		: base_t(arg1), m_op(op), m_arg2v(arg2v)
		{ }

		LMAT_ENSURE_INLINE const Op& op() const
		{
			return m_op;
		}

		LMAT_ENSURE_INLINE const T2& arg2_value() const
		{
			return m_arg2v;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		Op m_op;
		const T2 m_arg2v;
	};



	/********************************************
	 *
	 *  Expression mapping and construction
	 *
	 ********************************************/

	// spec classes

	template<typename Op>
	struct ewise_t
	{
		Op op;

		LMAT_ENSURE_INLINE
		ewise_t(const Op& o) : op(o) { }
	};

	template<typename Op, typename T>
	struct ewise_fix1st_t
	{
		Op op;
		const T& val;

		LMAT_ENSURE_INLINE
		ewise_fix1st_t(const Op& o, const T& v) : op(o), val(v) { }
	};


	template<typename Op, typename T>
	struct ewise_fix2nd_t
	{
		Op op;
		const T& val;

		LMAT_ENSURE_INLINE
		ewise_fix2nd_t(const Op& o, const T& v) : op(o), val(v) { }
	};


	// spec construction functions

	template<typename Op>
	LMAT_ENSURE_INLINE
	inline ewise_t<Op> ewise(const Op& op)
	{
		return ewise_t<Op>(op);
	}

	template<typename Op, typename T>
	LMAT_ENSURE_INLINE
	inline ewise_fix1st_t<Op, T> ewise_fix1st(const Op& op, const T& v)
	{
		return ewise_fix1st_t<Op, T>(op, v);
	}

	template<typename Op, typename T>
	LMAT_ENSURE_INLINE
	inline ewise_fix2nd_t<Op, T> ewise_fix2nd(const Op& op, const T& v)
	{
		return ewise_fix2nd_t<Op, T>(op, v);
	}


	// verifiers

	template<typename Op, class Arg>
	struct unary_expr_verifier<ewise_t<Op>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<typename Op, class Arg1, class Arg2>
	struct binary_expr_verifier<ewise_t<Op>, Arg1, Arg2>
	{
		static const bool value = is_mat_xpr<Arg1>::value && is_mat_xpr<Arg2>::value;
	};

	template<typename Op, typename T, class Arg>
	struct unary_expr_verifier<ewise_fix1st_t<Op, T>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<typename Op, typename T, class Arg>
	struct unary_expr_verifier<ewise_fix2nd_t<Op, T>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	// maps

	template<typename Op, typename Arg_HP, class Arg>
	struct unary_expr_map<ewise_t<Op>, Arg_HP, Arg>
	{
		typedef unary_ewise_expr<Op, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Op>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.op, arg_fwd);
		}
	};

	template<typename Op, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct binary_expr_map<ewise_t<Op>, Arg1_HP, Arg1, Arg2_HP, Arg2>
	{
		typedef binary_ewise_expr<Op, Arg1_HP, Arg1, Arg2_HP, Arg2> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Op>& spec,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		{
			return type(spec.op, arg1_fwd, arg2_fwd);
		}
	};

	template<typename Op, typename T1, typename Arg2_HP, class Arg2>
	struct unary_expr_map<ewise_fix1st_t<Op, T1>, Arg2_HP, Arg2>
	{
		typedef binary_fix1st_ewise_expr<Op, T1, Arg2_HP, Arg2> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_fix1st_t<Op, T1>& spec,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		{
			return type(spec.op, spec.val, arg2_fwd);
		}
	};

	template<typename Op, typename Arg1_HP, class Arg1, typename T2>
	struct unary_expr_map<ewise_fix2nd_t<Op, T2>, Arg1_HP, Arg1>
	{
		typedef binary_fix2nd_ewise_expr<Op, Arg1_HP, Arg1, T2> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_fix2nd_t<Op, T2>& spec,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd)
		{
			return type(spec.op, arg1_fwd, spec.val);
		}
	};


	// ewise expression construction

	template<typename Op, typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_t<Op>,
		ref_arg_t, Arg
	>::type
	ewise(const Op& op, const IMatrixXpr<Arg, T>& arg)
	{
		return make_expr( ewise(op), ref_arg(arg.derived()) );
	}

	template<typename Op, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Op>,
		ref_arg_t, Arg1,
		ref_arg_t, Arg2
	>::type
	ewise(  const Op& op,
			const IMatrixXpr<Arg1, T1>& arg1,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr(ewise(op),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}

	template<typename Op, typename T1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_fix1st_t<Op, T1>,
		ref_arg_t, Arg2
	>::type
	ewise_fix1st(  const Op& op,
			const T1& arg1v,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr( ewise_fix1st(op, arg1v), ref_arg(arg2.derived()) );
	}


	template<typename Op, typename T1, class Arg1, typename T2>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_fix2nd_t<Op, T2>,
		ref_arg_t, Arg1
	>::type
	ewise_fix2nd(  const Op& op,
			const IMatrixXpr<Arg1, T1>& arg1,
			const T2& arg2v )
	{
		return make_expr( ewise_fix2nd(op, arg2v), ref_arg(arg1.derived()) );
	}

}


/********************************************
 *
 *  Useful macros
 *
 ********************************************/

#define LMAT_DEFINE_UNARY_MATFUNCTION( Fname, Op ) \
		template<typename T, class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_t<Op>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, T>& arg) \
		{ return ewise( Op(), arg ); }

#define LMAT_DEFINE_BINARY_MATFUNCTION( Fname, Op ) \
		template<typename T, class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Op>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise( Op(), arg1, arg2 ); } \
		template<typename T, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_fix1st_t<Op, T>, ref_arg_t, Arg2>::type \
		Fname (const T& arg1v, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise_fix1st( Op(), arg1v, arg2 ); } \
		template<typename T, class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_fix2nd_t<Op, T>, ref_arg_t, Arg1>::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const T& arg2v) \
		{ return ewise_fix2nd( Op(), arg1, arg2v ); }



#endif 




