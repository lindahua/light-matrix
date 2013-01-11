/**
 * @file special.h
 *
 * @brief Special functions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_SPECIAL_H_
#define LIGHTMAT_SPECIAL_H_

#include <light_mat/math/math_base.h>
#include <cmath>

#include "internal/norminv_impl.h"

namespace lmat { namespace math {

#ifdef LMAT_HAS_CXX11_MATH

	using std::lgamma;
	using std::tgamma;

	using std::erf;
	using std::erfc;

#endif

	// norminv

	inline float norminv(float x)
	{
		return internal::norminv_impl<float>::eval(x);
	}

	inline double norminv(double x)
	{
		return internal::norminv_impl<double>::eval(x);
	}


} }

#endif
