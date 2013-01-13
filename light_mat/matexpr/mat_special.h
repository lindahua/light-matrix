/**
 * @file mat_special.h
 *
 * @brief Special functions on matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_SPECIAL_H_
#define LIGHTMAT_MAT_SPECIAL_H_

#include <light_mat/matexpr/matfun_base.h>
#include <light_mat/math/special_functors.h>

namespace lmat
{
	// gauss related

	_LMAT_DEFINE_RMATFUN( erf, 1 )
	_LMAT_DEFINE_RMATFUN( erfc, 1 )
	_LMAT_DEFINE_RMATFUN( norminv, 1 )

	// gamma related

	_LMAT_DEFINE_RMATFUN( lgamma, 1 )
	_LMAT_DEFINE_RMATFUN( tgamma, 1 )

}

#endif
