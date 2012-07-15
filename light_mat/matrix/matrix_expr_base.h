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

	template<class Expr, typename T>
	class embed_mat
	{
	public:
		LMAT_ENSURE_INLINE
		explicit embed_mat(const IMatrixXpr<Expr, T>& expr)
		: m_expr(expr.derived())
		{
		}

		LMAT_ENSURE_INLINE
		const Expr& get() const
		{
			return m_expr;
		}

	private:
		const Expr& m_expr;
	};


	template<class Expr, typename T>
	LMAT_ENSURE_INLINE
	embed_mat<Expr, T> embed(const IMatrixXpr<Expr, T>& expr)
	{
		return embed_mat<Expr, T>(expr);
	}
}


// Useful macros

#define LMAT_DECLARE_UNARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class Arg, bool EmbedArg> \
	struct map_prefix##_expr_map { \
		typedef typename unary_ewise_expr_map< \
				functor_name<typename matrix_traits<Arg>::value_type>, Arg, EmbedArg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP( map_prefix, functor_name ) \
	template<class LArg, class RArg, bool EmbedLArg, bool EmbedRArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg, EmbedLArg, EmbedRArg>::type type; };

#define LMAT_DECLARE_BINARY_TYPE_MAP_EX( map_prefix, functor_name ) \
	template<class LArg, class RArg, bool EmbedLArg, bool EmbedRArg> \
	struct map_prefix##_expr_map { \
		typedef typename binary_ewise_expr_map< \
				functor_name<typename binary_value_type<LArg, RArg>::type>, LArg, RArg, EmbedLArg, EmbedRArg>::type type; }; \
	template<class RArg, bool EmbedRArg> \
	struct map_prefix##_fix1_expr_map { \
		typedef typename binary_fix1_ewise_expr_map< \
				functor_name<typename matrix_traits<RArg>::value_type>, RArg, EmbedRArg>::type type; }; \
	template<class LArg, bool EmbedLArg> \
	struct map_prefix##_fix2_expr_map { \
		typedef typename binary_fix2_ewise_expr_map< \
				functor_name<typename matrix_traits<LArg>::value_type>, LArg, EmbedLArg>::type type; };


#define LMAT_DEFINE_UNARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class Arg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<Arg, false>::type \
	matfun_name(const IMatrixXpr<Arg, T>& A) { \
		return ewise(functor_name<T>(), A.derived()); } \
	template<typename T, class Arg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<Arg, true>::type \
	matfun_name(const embed_mat<Arg, T>& A) { \
	return ewise(functor_name<T>(), A); }


#define LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg, false, false>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), A.derived(), B.derived()); } \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg, true, false>::type \
	matfun_name(const embed_mat<LArg, T>& A, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), A, B.derived()); } \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg, false, true>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const embed_mat<RArg, T>& B) { \
		return ewise(functor_name<T>(), A.derived(), B); } \
	template<typename T, class LArg, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_expr_map<LArg, RArg, true, true>::type \
	matfun_name(const embed_mat<LArg, T>& A, const embed_mat<RArg, T>& B) { \
		return ewise(functor_name<T>(), A, B); }


#define LMAT_DEFINE_BINARY_MATFUNCTION_EX( map_prefix, matfun_name, functor_name ) \
	LMAT_DEFINE_BINARY_MATFUNCTION( map_prefix, matfun_name, functor_name ) \
	template<typename T, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix1_expr_map<RArg, false>::type \
	matfun_name(const T& a, const IMatrixXpr<RArg, T>& B) { \
		return ewise(functor_name<T>(), a, B.derived()); } \
	template<typename T, class RArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix1_expr_map<RArg, true>::type \
	matfun_name(const T& a, const embed_mat<RArg, T>& B) { \
		return ewise(functor_name<T>(), a, B); } \
	template<typename T, class LArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix2_expr_map<LArg, false>::type \
	matfun_name(const IMatrixXpr<LArg, T>& A, const T& b) { \
		return ewise(functor_name<T>(), A.derived(), b); } \
	template<typename T, class LArg> \
	LMAT_ENSURE_INLINE \
	inline typename map_prefix##_fix2_expr_map<LArg, true>::type \
	matfun_name(const embed_mat<LArg, T>& A, const T& b) { \
		return ewise(functor_name<T>(), A, b); }

#endif






