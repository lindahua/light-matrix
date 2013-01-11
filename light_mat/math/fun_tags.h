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

// macros to declare fun tags


namespace lmat { namespace ftags {


	// arithmetic

	struct add_ { };
	struct sub_ { };
	struct mul_ { };
	struct div_ { };
	struct neg_ { };

	struct abs_ { };
	struct sqr_ { };
	struct cube_ { };

	struct max_ { };
	struct min_ { };
	struct clamp_ { };

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


	// real math

	struct fma_ { };
	struct rcp_ { };
	struct sqrt_ { };
	struct rsqrt_ { };
	struct pow_ { };

	struct floor_ { };
	struct ceil_ { };

	struct exp_ { };
	struct log_ { };
	struct log10_ { };
	struct xlogx_ { };
	struct xlogy_ { };

	struct sin_ { };
	struct cos_ { };
	struct tan_ { };

	struct asin_ { };
	struct acos_ { };
	struct atan_ { };
	struct atan2_ { };

	struct sinh_ { };
	struct cosh_ { };
	struct tanh_ { };

	// C++11 real math

	struct cbrt_ { };
	struct hypot_ { };

	struct round_ { };
	struct trunc_ { };

	struct exp2_ { };
	struct log2_ { };
	struct expm1_ { };
	struct log1p_ { };

	struct asinh_ { };
	struct acosh_ { };
	struct atanh_ { };

	// numeric predicates

	struct signbit_ { };
	struct isfinite_ { };
	struct isinf_ { };
	struct isnan_ { };

	struct cond_ { };

	// special functions

	struct erf_ { };
	struct erfc_ { };
	struct norminv_ { };

	struct lgamma_ { };
	struct tgamma_ { };
	struct psi_ { };

} }

#endif 
