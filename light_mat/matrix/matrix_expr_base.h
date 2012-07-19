/**
 * @file matrix_expr_base.h
 *
 * The basis of all matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_EXPR_BASE_H_
#define LIGHTMAT_MATRIX_EXPR_BASE_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{
	template<class Expr>
	struct matrix_traits<embed_mat<Expr> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Expr>::value;
		static const int compile_time_num_cols = ct_cols<Expr>::value;

		static const bool is_readonly = true;
		typedef typename matrix_traits<Expr>::value_type value_type;
	};

	template<class Expr>
	class embed_mat : public IMatrixXpr<embed_mat<Expr>, typename matrix_traits<Expr>::value_type>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit embed_mat(const Expr& expr)
		: m_expr(expr.derived())
		{
		}

		LMAT_ENSURE_INLINE
		const Expr& internal() const
		{
			return m_expr;
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_expr.nelems();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return m_expr.size();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_expr.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_expr.ncolumns();
		}

	private:
		const Expr& m_expr;
	};


	template<class Expr, typename T>
	LMAT_ENSURE_INLINE
	embed_mat<Expr> embed(const IMatrixXpr<Expr, T>& expr)
	{
		return embed_mat<Expr>(expr.derived());
	}


	template<class Expr>
	class obj_wrapper<embed_mat<Expr> >
	{
	public:
		LMAT_ENSURE_INLINE
		obj_wrapper(const embed_mat<Expr>& e)
		: m_expr(e.internal())
		{
		}

		LMAT_ENSURE_INLINE
		const Expr& get() const
		{
			return m_expr;
		}
	private:
		Expr m_expr;
	};


	template<class Expr>
	struct unwrapped_expr
	{
		typedef Expr type;
	};

	template<class Expr>
	struct unwrapped_expr<embed_mat<Expr> >
	{
		typedef Expr type;
	};

}


// Useful macros

#define LMAT_DECLARE_UNARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class Arg> \
	struct map_prefix##_expr_map { \
		typedef typename unary_ewise_expr_map< \
				functor_name<typename matrix_traits<Arg>::value_type>, Arg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class LArg, class RArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP_EX( map_prefix, functor_name ) \
	template<class LArg, class RArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg>::type type; }; \
	template<class RArg> \
	struct map_prefix##_fix1_expr_map { \
		typedef typename binary_fix1_ewise_expr_map< \
				functor_name<typename matrix_traits<RArg>::value_type>, RArg>::type type; }; \
	template<class LArg> \
	struct map_prefix##_fix2_expr_map { \
		typedef typename binary_fix2_ewise_expr_map< \
				functor_name<typename matrix_traits<LArg>::value_type>, LArg>::type type; };


#define LMAT_DEFINE_UNARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class Arg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<Arg>::type \
	matfun_name(const IMatrixXpr<Arg, T>& A) { \
		return ewise(functor_name<T>(), A.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), A.derived(), B.derived()); }


#define LMAT_DEFINE_BINARY_MATFUNCTION_EX( map_prefix, matfun_name, functor_name ) \
	LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix1_expr_map<RArg>::type \
	matfun_name(const T& a, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), a, B.derived()); } \
	template<typename T, class LArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix2_expr_map<LArg>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const T& b) { \
		return ewise(functor_name<T>(), A.derived(), b); }

#endif






