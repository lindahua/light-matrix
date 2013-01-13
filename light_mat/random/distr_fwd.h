/**
 * @file prng_fwd.h
 *
 * @brief Forward header for all distributions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DISTR_FWD_H_
#define LIGHTMAT_DISTR_FWD_H_

#include <light_mat/random/sfmt.h>
#include <light_mat/math/math_base.h>
#include <light_mat/math/functor_base.h>

namespace lmat { namespace random {


	// tags to indicate PRNG methods

	struct naive_ { };
	struct basic_ { };
	struct icdf_ { };
	struct box_muller_ { };
	struct marsaglia_ { };
	struct ziggurat_ { };
	struct huffman_ { };

	// discrete distributions

	template<typename TI=uint32_t> class std_uniform_int_distr;
	template<typename TI=uint32_t> class uniform_int_distr;

	class std_bernoulli_distr;
	class bernoulli_distr;

	template<typename TI=uint32_t, typename Method=naive_> class binomial_distr;
	template<typename TI=uint32_t, typename Method=naive_> class geometric_distr;
	template<typename TI=uint32_t, typename Method=naive_> class negative_binomial_distr;
	template<typename TI=uint32_t, typename Method=naive_> class poisson_distr;

	template<typename TI=uint32_t, typename Method=naive_> class discrete_distr;

	// real-value distributions

	template<typename T=double> class std_uniform_real_distr;
	template<typename T=double> class uniform_real_distr;

	template<typename T=double> class std_exponential_distr;
	template<typename T=double> class exponential_distr;

	template<typename T=double, typename Method=icdf_> class std_normal_distr;
	template<typename T=double, typename Method=icdf_> class normal_distr;

	template<typename T=double, typename Method=basic_> class std_gamma_distr;
	template<typename T=double, typename Method=basic_> class gamma_distr;


} }

#endif
