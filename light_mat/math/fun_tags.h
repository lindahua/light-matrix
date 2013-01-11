/**
 * @file fun_tags.h
 *
 * Function tags
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUN_TAGS_H_
#define LIGHTMAT_FUN_TAGS_H_

#include <light_mat/common/prim_types.h>
#include <light_mat/common/mask_type.h>

namespace lmat { namespace ftags {

	// arithmetic

	struct add_ { };
	struct sub_ { };
	struct mul_ { };
	struct div_ { };
	struct neg_ { };
	struct fma_ { };

	struct max_ { };
	struct min_ { };
	struct clamp_ { };
	struct cond_ { };

	// comparison

	struct eq_ { };
	struct ne_ { };
	struct ge_ { };
	struct gt_ { };
	struct le_ { };
	struct lt_ { };

	// logical

	struct logical_not_ { };
	struct logical_and_ { };
	struct logical_or_ { };
	struct logical_eq_ { };
	struct logical_ne_ { };

	// simple power functions

	struct abs_ { };
	struct sqr_ { };
	struct cube_ { };

	struct rcp_ { };
	struct sqrt_ { };
	struct rsqrt_ { };

	// rounding

	struct floor_ { };
	struct ceil_ { };
	struct round_ { };
	struct trunc_ { };

	// numeric predicates

	struct signbit_ { };
	struct isfinite_ { };
	struct isinf_ { };
	struct isnan_ { };

	// ******************************

	// power functions

	struct pow_ { };
	struct cbrt_ { };
	struct hypot_ { };

	// exp & log

	struct exp_ { };
	struct log_ { };
	struct log10_ { };
	struct xlogx_ { };
	struct xlogy_ { };

	struct exp2_ { };
	struct log2_ { };
	struct exp10_ { };
	struct expm1_ { };
	struct log1p_ { };

	// trigonometry

	struct sin_ { };
	struct cos_ { };
	struct tan_ { };

	struct asin_ { };
	struct acos_ { };
	struct atan_ { };
	struct atan2_ { };

	// hyperbolic

	struct sinh_ { };
	struct cosh_ { };
	struct tanh_ { };

	struct asinh_ { };
	struct acosh_ { };
	struct atanh_ { };

	// ******************************

	// special functions

	struct erf_ { };
	struct erfc_ { };
	struct erfinv_ { };
	struct erfcinv_ { };
	struct norminv_ { };

	struct lgamma_ { };
	struct tgamma_ { };
	struct psi_ { };

} }

#endif 
