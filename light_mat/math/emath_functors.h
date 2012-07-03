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
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( min_fun, (math::min), true )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( max_fun, (math::max), true )

	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( pow_fun,  math::pow, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( floor_fun, math::floor, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( ceil_fun,  math::ceil, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( exp_fun,   math::exp, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log_fun,   math::log, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log10_fun, math::log10, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( sin_fun, math::sin, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cos_fun, math::cos, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( tan_fun, math::tan, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( asin_fun, math::asin, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( acos_fun, math::acos, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( atan_fun, math::atan, true )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( atan2_fun, math::atan2, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( sinh_fun, math::sinh, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cosh_fun, math::cosh, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( tanh_fun, math::tanh, true )

#ifdef LMAT_HAS_CXX11_MATH

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( round_fun, math::round, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( trunc_fun, math::trunc, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( cbrt_fun,  math::cbrt, true )
	LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( hypot_fun, math::hypot, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( exp2_fun, math::exp2, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log2_fun, math::log2, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( expm1_fun, math::expm1, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( log1p_fun, math::log1p, true )

	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( asinh_fun, math::asinh, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( acosh_fun, math::acosh, true )
	LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( atanh_fun, math::atanh, true )

#endif

}

#endif /* EMATH_FUNCTORS_H_ */
