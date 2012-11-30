/**
 * @file matrix_cast.h
 *
 * Value-type casting for matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_CAST_H_
#define LIGHTMAT_MATRIX_CAST_H_

#include <light_mat/matexpr/matrix_ewise_eval.h>

#define LMAT_DEFINE_CASTING_MATFUNCTION( fname, ttype ) \
	template<class SExpr, typename S> \
	LMAT_ENSURE_INLINE \
	inline typename unary_expr_map<ewise_t<scast_t<ttype> >, ref_arg_t, SExpr >::type \
	fname(const IMatrixXpr<SExpr, S>& sexpr) \
	{ return cast(sexpr, type<ttype>()); }


namespace lmat
{
	// functor definitions

	template<typename T>
	struct scast_t { };

	template<typename T>
	struct is_unary_op<scast_t<T> >
	{
		static const bool value = true;
	};

	template<typename S, typename T>
	struct unary_op_result<scast_t<T>, S>
	{
		typedef T type;
	};

	template<typename S, typename T>
	struct scast_fun
	{
		typedef S source_type;
		typedef T target_type;

		LMAT_ENSURE_INLINE
		scast_fun() { }

		LMAT_ENSURE_INLINE
		scast_fun( scast_t<T> ) { }

		LMAT_ENSURE_INLINE
		T operator() (const S& s) const
		{
			return static_cast<T>(s);
		}
	};

	template<typename S, typename T>
	struct unary_op_fun<scast_t<T>, scalar_ker, S>
	{
		typedef scast_fun<S, T> type;
	};


	// matrix functions

	template<class SExpr, typename S, typename T>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<
		ewise_t<scast_t<T> >,
		ref_arg_t, SExpr >::type
	cast(const IMatrixXpr<SExpr, S>& sexpr, type<T> )
	{
		return ewise(scast_t<T>(), sexpr.derived());
	}

	LMAT_DEFINE_CASTING_MATFUNCTION( to_bool, bool )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_f32, float )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_f64, double )

	LMAT_DEFINE_CASTING_MATFUNCTION( to_i32, int32_t )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_u32, uint32_t )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_i16, int16_t )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_u16, uint16_t )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_i8,  int8_t )
	LMAT_DEFINE_CASTING_MATFUNCTION( to_u8,  uint8_t )

}

#endif /* MATRIX_CAST_H_ */
