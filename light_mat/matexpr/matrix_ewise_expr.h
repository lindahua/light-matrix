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
#include <light_mat/matexpr/matexpr_tools.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expression class traits
	 *
	 ********************************************/

	template<typename Tag, class QLst>
	struct matrix_traits<ewise_expr<Tag, QLst> >
	{
		typedef typename meta::argtype_list_<QLst>::type arg_list;

		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::common_nrows<arg_list>::value;
		static const int ct_num_cols = meta::common_ncols<arg_list>::value;

		static const bool is_readonly = true;

		typedef typename meta::transform_<meta::value_type_of, arg_list>::type arg_value_types;
		typedef typename op_result<Tag, arg_value_types>::type value_type;
		typedef typename meta::common_domain<arg_list>::type domain;
	};


	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<typename Tag, class QLst>
	class ewise_expr
	: public IMatrixXpr<ewise_expr<Tag, QLst>,
		typename matrix_traits< ewise_expr<Tag, QLst> >::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_op_tag<Tag>::value, "Tag must be an operation tag");
#endif

		typedef typename meta::argtype_list_<QLst>::type arg_list;

		typedef matrix_shape<
				meta::common_nrows<arg_list>::value,
				meta::common_ncols<arg_list>::value > shape_type;

	public:
		LMAT_ENSURE_INLINE
		ewise_expr(const Tag& t, const tied_forwarder<QLst>& arg_fwd)
		: m_tag(t)
		, m_shape( common_nrows_tfwd(arg_fwd), common_ncols_tfwd(arg_fwd) )
		, m_args(arg_fwd) { }

		LMAT_ENSURE_INLINE const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE const qarg_list<QLst>& args() const
		{
			return m_args;
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
		qarg_list<QLst> m_args;
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

	template<typename Tag, class ArgList>
	struct expr_verifier<ewise_t<Tag>, ArgList>
	{
		static const bool value = meta::have_compatible_shapes<ArgList>::value;
	};

	// maps

	template<typename Tag, class QLst>
	struct expr_map<ewise_t<Tag>, QLst>
	{
		typedef ewise_expr<Tag, QLst> type;

		LMAT_ENSURE_INLINE
		static type get(const ewise_t<Tag>& spec,
				const tied_forwarder<QLst>& tfwd)
		{
			return type(spec.tag, tfwd);
		}
	};


	// ewise expression construction

	template<typename Tag, typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename expr_map< ewise_t<Tag>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type
	ewise(const Tag& tag, const IMatrixXpr<Arg, T>& arg)
	{
		typedef LMAT_TYPELIST_1(CRefArg<Arg>) QLst;
		return make_expr( ewise(tag), tied_forwarder<QLst>(arg.derived()) );
	}

	template<typename Tag, typename T1, class Arg1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename expr_map< ewise_t<Tag>, LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2>) >::type
	ewise(  const Tag& tag,
			const IMatrixXpr<Arg1, T1>& arg1,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		typedef LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2>) QLst;
		return make_expr( ewise(tag), tied_forwarder<QLst>(arg1.derived(), arg2.derived()) );
	}

	template<typename Tag, typename T1, typename T2, class Arg2>
	LMAT_ENSURE_INLINE
	inline typename expr_map< ewise_t<Tag>,
		LMAT_TYPELIST_2( CpyArg<scalar_expr<T1> >, CRefArg<Arg2> ) >::type
	ewise(  const Tag& tag,
			const scalar_expr<T1>& arg1v,
			const IMatrixXpr<Arg2, T2>& arg2 )
	{
		typedef LMAT_TYPELIST_2( CpyArg<scalar_expr<T1> >, CRefArg<Arg2> ) QLst;
		return make_expr( ewise(tag), tied_forwarder<QLst>(arg1v, arg2.derived()) );
	}


	template<typename Tag, typename T1, class Arg1, typename T2>
	LMAT_ENSURE_INLINE
	inline typename expr_map<ewise_t<Tag>,
		LMAT_TYPELIST_2( CRefArg<Arg1>, CpyArg<scalar_expr<T2> > ) >::type
	ewise(  const Tag& tag,
			const IMatrixXpr<Arg1, T1>& arg1,
			const scalar_expr<T2>& arg2v )
	{
		typedef LMAT_TYPELIST_2( CRefArg<Arg1>, CpyArg<scalar_expr<T2> > ) QLst;
		return make_expr( ewise(tag), tied_forwarder<QLst>(arg1.derived(), arg2v) );
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
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type \
		Fname (const IMatrixXpr<Arg, T>& arg) \
		{ return ewise( Tag(), arg ); }

#define LMAT_DEFINE_BINARY_MATFUNCTION( Fname, Tag ) \
		template<typename T, class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2> ) >::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise( Tag(), arg1, arg2 ); } \
		template<typename T, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CpyArg<scalar_expr<T> >, CRefArg<Arg2> ) >::type \
		Fname (const T& arg1v, const IMatrixXpr<Arg2, T>& arg2) \
		{ return ewise( Tag(), make_scalar(arg1v), arg2 ); } \
		template<typename T, class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CRefArg<Arg1>, CpyArg<scalar_expr<T> > ) >::type \
		Fname (const IMatrixXpr<Arg1, T>& arg1, const T& arg2v) \
		{ return ewise( Tag(), arg1, make_scalar(arg2v) ); }

// of specific input value type

#define LMAT_DEFINE_UNARY_MATFUNCTION_S( Fname, Tag, SType ) \
		template<class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type \
		Fname (const IMatrixXpr<Arg, SType>& arg) \
		{ return ewise( Tag(), arg ); }

#define LMAT_DEFINE_BINARY_MATFUNCTION_S( Fname, Tag, SType ) \
		template<class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CRefArg<Arg1>, CRefArg<Arg2> ) >::type \
		Fname (const IMatrixXpr<Arg1, SType>& arg1, const IMatrixXpr<Arg2, SType>& arg2) \
		{ return ewise( Tag(), arg1, arg2 ); } \
		template<class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CpyArg<scalar_expr<SType> >, CRefArg<Arg2> ) >::type \
		Fname (const SType& arg1v, const IMatrixXpr<Arg2, SType>& arg2) \
		{ return ewise( Tag(), make_scalar(arg1v), arg2 ); } \
		template<class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename expr_map<ewise_t<Tag>, LMAT_TYPELIST_2(CRefArg<Arg1>, CpyArg<scalar_expr<SType> > ) >::type \
		Fname (const IMatrixXpr<Arg1, SType>& arg1, const SType& arg2v) \
		{ return ewise( Tag(), arg1, make_scalar(arg2v) ); }

#endif 




