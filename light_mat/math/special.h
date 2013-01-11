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

namespace lmat { namespace math {

#ifdef LMAT_HAS_CXX11_MATH

	using std::lgamma;
	using std::tgamma;

	using std::erf;
	using std::erfc;

#endif

} }

#endif
