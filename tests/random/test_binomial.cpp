/**
 * @file test_binomial.cpp
 *
 * @brief Unit testing of binomial distributions
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/binomial_distr.h>

default_rand_stream rstream;
const index_t N = 200000;

SIMPLE_CASE( test_binomial_naive )
{
	uint32_t t = 5;
	const double p = 0.4;
	binomial_distr<uint32_t, naive_> distr(t, p);

	ASSERT_EQ( distr.t(), t );
	ASSERT_EQ( distr.p(), p );
	ASSERT_EQ( distr.mean(), double(t) * p );
	ASSERT_APPROX( distr.var(), double(t) * p * (1-p), 1.0e-15 );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, (index_t)t+1, ptol );
}

AUTO_TPACK( test_binomial )
{
	ADD_SIMPLE_CASE( test_binomial_naive )
}

