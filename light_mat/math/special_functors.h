/**
 * @file special_functors.h
 *
 * @brief Functors for special math functions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SPECIAL_FUNCTORS_H_
#define LIGHTMAT_SPECIAL_FUNCTORS_H_

#include <light_mat/math/basic_functors.h>
#include <light_mat/math/math_special.h>
#include <light_mat/math/simd_math.h>

namespace lmat
{
	// gauss related functions

	_LMAT_DEFINE_REAL_MATH_FUN( erf, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( erfc, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( norminv, 1 )

	// gamma related functions

	_LMAT_DEFINE_REAL_MATH_FUN( lgamma, 1 )
	_LMAT_DEFINE_REAL_MATH_FUN( tgamma, 1 )
}


#endif
