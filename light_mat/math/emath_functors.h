/**
 * @file emath_functors.h
 *
 * Elementary math functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EMATH_FUNCTORS_H_
#define LIGHTMAT_EMATH_FUNCTORS_H_

#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>

namespace lmat
{
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( min_fun, (math::min) )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( max_fun, (math::max) )

	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( pow_fun,  math::pow )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( floor_fun, math::floor )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( ceil_fun,  math::ceil )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( exp_fun,   math::exp )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log_fun,   math::log )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log10_fun, math::log10 )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( sin_fun, math::sin )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cos_fun, math::cos )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( tan_fun, math::tan )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( asin_fun, math::asin )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( acos_fun, math::acos )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( atan_fun, math::atan )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( atan2_fun, math::atan2 )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( sinh_fun, math::sinh )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cosh_fun, math::cosh )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( tanh_fun, math::tanh )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( round_fun, math::round )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( trunc_fun, math::trunc )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cbrt_fun,  math::cbrt )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( hypot_fun, math::hypot )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( exp2_fun, math::exp2 )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log2_fun, math::log2 )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( expm1_fun, math::expm1 )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log1p_fun, math::log1p )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( asinh_fun, math::asinh )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( acosh_fun, math::acosh )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( atanh_fun, math::atanh )

#endif

}

#endif /* EMATH_FUNCTORS_H_ */
