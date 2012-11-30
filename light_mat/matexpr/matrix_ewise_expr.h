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

	template<typename Tag, typename Arg_HP, class Arg>
	struct matrix_traits<unary_ewise_expr<Tag, Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value;
		static const int ct_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type arg_value_type;

		typedef typename unary_op_result<Tag, arg_value_type>::type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Tag, typename Arg1_HP, class Arg1, class Arg2_HP, class Arg2>
	struct matrix_traits<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = common_ctrows<Arg1, Arg2>::value;
		static const int ct_num_cols = common_ctcols<Arg1, Arg2>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg1>::value_type arg1_value_type;
		typedef typename matrix_traits<Arg2>::value_type arg2_value_type;

		typedef typename binary_op_result<Tag, arg1_value_type, arg2_value_type>::type value_type;
		typedef typename binary_domain<Arg1, Arg2>::type domain;
	};

	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<typename Tag, typename Arg_HP, class Arg>
	class unary_ewise_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<unary_ewise_expr<Tag, Arg_HP, Arg>,
		typename matrix_traits<unary_ewise_expr<Tag, Arg_HP, Arg> >::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_unary_op<Tag>::value, "Tag must be a unary operation");
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;

		LMAT_ENSURE_INLINE
		unary_ewise_expr(const Tag& t,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd), m_tag(t) { }

		LMAT_ENSURE_INLINE const Tag& tag() const
		{
			return m_tag;
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
		Tag m_tag;
	};


	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	class binary_ewise_expr
	: public binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2>
	, public IMatrixXpr<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2>,
		typename matrix_traits<binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2> >::value_type >
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_binary_op<Tag>::value, "Tag must be a binary operation.");
#endif

	public:
		typedef binary_expr_base<Arg1_HP, Arg1, Arg2_HP, Arg2> base_t;
		typedef typename common_shape<Arg1, Arg2>::type shape_type;

		LMAT_ENSURE_INLINE
		binary_ewise_expr(const Tag& t,
				const arg_forwarder<Arg1_HP, Arg1>& arg1,
				const arg_forwarder<Arg2_HP, Arg2>& arg2)
		: base_t(arg1, arg2)
		, m_shape(get_common_nrows(arg1.arg, arg2.arg), get_common_ncolumns(arg1.arg, arg2.arg))
		, m_tag(t)
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
		shape_type m_shape;
		Tag m_tag;
	};


	/********************************************
	 *
	 *  Expression mapping and construction
	 *
	 ********************************************/

	// spec classes

	template<typename Tag>
	struct ewise_t
	{
		Tag tag;

		LMAT_ENSURE_INLINE
		ewise_t(const Tag& t) : tag(t) { }
	};

	template<typename Tag>
	LMAT_ENSURE_INLINE
	inline ewise_t<Tag> ewise(const Tag& op)
	{
		return ewise_t<Tag>(op);
	}

	// verifiers

	template<typename Tag, class Arg>
	struct unary_expr_verifier<ewise_t<Tag>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<typename Tag, class Arg1, class Arg2>
	struct binary_expr_verifier<ewise_t<Tag>, Arg1, Arg2>
	{
		static const bool value =
				(is_mat_xpr<Arg1>::value || is_mat_xpr<Arg2>::value) &&
				has_compatible_ct_size<Arg1, Arg2>::value;
	};

	// maps

	template<typename Tag, typename Arg_HP, class Arg>
	struct unary_expr_map<ewise_t<Tag>, Arg_HP, Arg>
	{
		typedef unary_ewise_expr<Tag, Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Tag>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(spec.tag, arg_fwd);
		}
	};

	template<typename Tag, typename Arg1_HP, class Arg1, typename Arg2_HP, class Arg2>
	struct binary_expr_map<ewise_t<Tag>, Arg1_HP, Arg1, Arg2_HP, Arg2>
	{
		typedef binary_ewise_expr<Tag, Arg1_HP, Arg1, Arg2_HP, Arg2> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Tag>& spec,
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
		{
			return type(spec.tag, arg1_fwd, arg2_fwd);
		}
	};



	// ewise expression construction

	template<typename Tag, typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<ewise_t<Tag>,
		ref_arg_t, Arg
	>::type
	ewise(const Tag& tag, const IMatrixXpr<Arg, T>& arg)
	{
		return make_expr( ewise(tag), ref_arg(arg.derived()) );
	}

	template<typename Tag, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Tag>,
		ref_arg_t, Arg1,
		ref_arg_t, Arg2
	>::type
	ewise(  const Tag& tag,
			const IMatrixXpr<Arg1, T1>& arg1,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr(ewise(tag),
				ref_arg(arg1.derived()),
				ref_arg(arg2.derived()));
	}

	template<typename Tag, typename T1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Tag>,
		copy_arg_t, scalar_expr<T1>,
		ref_arg_t, Arg2
	>::type
	ewise(  const Tag& tag,
			const scalar_expr<T1>& arg1v,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		return make_expr( ewise(tag),
				copy_arg(arg1v),
				ref_arg(arg2.derived()) );
	}


	template<typename Tag, typename T1, class Arg1, typename T2>
	LMAT_ENSURE_INLINE
	inline typename binary_expr_map<ewise_t<Tag>,
		ref_arg_t, Arg1,
		copy_arg_t, scalar_expr<T2>
	>::type
	ewise(  const Tag& tag,
			const IMatrixXpr<Arg1, T1>& arg1,
			const scalar_expr<T2>& arg2v )
	{
		return make_expr( ewise(tag),
				ref_arg(arg1.derived()),
				copy_arg(arg2v) );
	}

}


/********************************************
 *
 *  Useful macros to define mat functions
 *
 ********************************************/

// generic input value type

#define LMAT_DEFINE_UNARY_MATFUNCTION( Fname, Tag ) \
		template<typename T, class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_t<Tag>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, T>& arg) \
		{ return ewise( Tag(), arg ); }

#define LMAT_DEFINE_BINARY_MATFUNCTION( Fname, Tag ) \
		template<typename T, class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise( Tag(), arg1, arg2 ); } \
		template<typename T, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, copy_arg_t, scalar_expr<T>, ref_arg_t, Arg2>::type \
		Fname (const T& arg1v, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise( Tag(), make_scalar(arg1v), arg2 ); } \
		template<typename T, class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, ref_arg_t, Arg1, copy_arg_t, scalar_expr<T> >::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const T& arg2v) \
		{ return ewise( Tag(), arg1, make_scalar(arg2v) ); }

// of specific input value type

#define LMAT_DEFINE_UNARY_MATFUNCTION_S( Fname, Tag, SType ) \
		template<class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_t<Tag>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, SType>& arg) \
		{ return ewise( Tag(), arg ); }

#define LMAT_DEFINE_BINARY_MATFUNCTION_S( Fname, Tag, SType ) \
		template<class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, SType>& arg1, const IMatrixXpr<Arg2, SType>& arg2) \
		{ return ewise( Tag(), arg1, arg2 ); } \
		template<class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, copy_arg_t, scalar_expr<SType>, ref_arg_t, Arg2>::type \
		Fname (const SType& arg1v, const IMatrixXpr<Arg2, SType>& arg2) \
		{ return ewise( Tag(), make_scalar(arg1v), arg2 ); } \
		template<class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Tag>, ref_arg_t, Arg1, copy_arg_t, scalar_expr<SType> >::type \
		Fname (const IMatrixXpr<Arg1, SType>& arg1, const SType& arg2v) \
		{ return ewise( Tag(), arg1, make_scalar(arg2v) ); }

#endif 




