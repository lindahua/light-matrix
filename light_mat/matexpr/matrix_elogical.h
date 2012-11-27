/*
 * @file matrix_elogical.h
 *
 * Element-wise logical operations
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ELOGICAL_H_
#define LIGHTMAT_MATRIX_ELOGICAL_H_

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include <light_mat/math/logical_functors.h>


#define LMAT_DEFINE_UNARY_LOGICAL_MATFUNCTION( Fname, Op ) \
		template<class Arg> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_t<Op>, ref_arg_t, Arg>::type \
		Fname (const IMatrixXpr<Arg, bool>& arg) \
		{ return ewise( Op(), arg ); }

#define LMAT_DEFINE_BINARY_LOGICAL_MATFUNCTION( Fname, Op ) \
		template<class Arg1, class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename binary_expr_map<ewise_t<Op>, ref_arg_t, Arg1, ref_arg_t, Arg2>::type \
		Fname (const IMatrixXpr<Arg1, bool>& arg1, const IMatrixXpr<Arg2, bool>& arg2) \
		{ return ewise( Op(), arg1, arg2 ); } \
		template<class Arg2> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_fix1st_t<Op, bool>, ref_arg_t, Arg2>::type \
		Fname (const bool& arg1v, const IMatrixXpr<Arg2, bool>& arg2) \
		{ return ewise_fix1st( Op(), arg1v, arg2 ); } \
		template<class Arg1> \
		LMAT_ENSURE_INLINE \
		inline typename unary_expr_map<ewise_fix2nd_t<Op, bool>, ref_arg_t, Arg1>::type \
		Fname (const IMatrixXpr<Arg1, bool>& arg1, const bool& arg2v) \
		{ return ewise_fix2nd( Op(), arg1, arg2v ); }


namespace lmat
{
	LMAT_DEFINE_UNARY_LOGICAL_MATFUNCTION ( operator ~, not_t )
	LMAT_DEFINE_BINARY_LOGICAL_MATFUNCTION( operator &, and_t )
	LMAT_DEFINE_BINARY_LOGICAL_MATFUNCTION( operator |, or_t )
	LMAT_DEFINE_BINARY_LOGICAL_MATFUNCTION( operator ^, xor_t )
}

#endif 
