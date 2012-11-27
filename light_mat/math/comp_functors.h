/*
 * @file comp_functors.h
 *
 * Element-wise comparison functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_COMP_FUNCTORS_H_
#define LIGHTMAT_COMP_FUNCTORS_H_

#include <light_mat/math/functor_base.h>


#define LMAT_DEFINE_COMPARISON_FUNCTOR( fname, opsym ) \
	template<typename T> struct fname##_fun { \
		LMAT_ENSURE_INLINE fname##_fun() { } \
		LMAT_ENSURE_INLINE fname##_fun( fname##_t ) { } \
		LMAT_ENSURE_INLINE bool operator() (const T& x, const T& y) const { return x opsym y; } \
	}; \
	LMAT_DEFINE_BINARY_PRED_FUNMAP( fname##_t, scalar_kernel_t, fname##_fun )

namespace lmat
{
	// declarations

	LMAT_DEFINE_BINARY_PRED_OP( eq_t )
	LMAT_DEFINE_BINARY_PRED_OP( ne_t )

	LMAT_DEFINE_BINARY_PRED_OP( gt_t )
	LMAT_DEFINE_BINARY_PRED_OP( ge_t )
	LMAT_DEFINE_BINARY_PRED_OP( lt_t )
	LMAT_DEFINE_BINARY_PRED_OP( le_t )

	/********************************************
	 *
	 *  scalar functors
	 *
	 ********************************************/

	LMAT_DEFINE_COMPARISON_FUNCTOR( eq, == )
	LMAT_DEFINE_COMPARISON_FUNCTOR( ne, != )
	LMAT_DEFINE_COMPARISON_FUNCTOR( ge, >= )
	LMAT_DEFINE_COMPARISON_FUNCTOR( gt, > )
	LMAT_DEFINE_COMPARISON_FUNCTOR( le, <= )
	LMAT_DEFINE_COMPARISON_FUNCTOR( lt, < )

}

#endif 
