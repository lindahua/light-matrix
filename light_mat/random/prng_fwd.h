/**
 * @file prng_fwd.h
 *
 * @brief Forward header for all Pseudo RNGs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PRNG_FWD_H_
#define LIGHTMAT_PRNG_FWD_H_

#include <light_mat/random/sfmt.h>

namespace lmat { namespace random {

	// forward declaration of distribution classes

	template<class Distr, typename Kind>
	struct distr_supp_simd
	{
		static const bool value = false;
	};

	// discrete distributions

	template<typename TI=uint32_t> class std_uniform_int_distr;
	template<typename TI=uint32_t> class uniform_int_distr;

	class std_bernoulli_distr;
	class bernoulli_distr;

	template<typename TI=uint32_t> class binomial_distr;
	template<typename TI=uint32_t> class negative_binomial_distr;
	template<typename TI=uint32_t> class geometric_distr;
	template<typename TI=uint32_t> class poisson_distr;

	template<typename TI=uint32_t> class discrete_distr;
	template<typename TI=uint32_t> class huffman_discrete_distr;

	// real-value distributions

	template<typename T=double> class std_uniform_real_distr;
	template<typename T=double> class uniform_real_distr;

	template<typename T=double> class std_normal_distr;
	template<typename T=double> class normal_distr;

	template<typename T=double> class std_exponential_distr;
	template<typename T=double> class exponential_distr;

	template<typename T=double> class std_gamma_distr;
	template<typename T=double> class gamma_distr;

	template<typename T=double> class std_lognormal_distr;
	template<typename T=double> class lognormal_distr;


} }

#endif
