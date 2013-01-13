/**
 * @file test_poissond.cpp
 *
 * @brief Unit testing of Poisson distribution
 *
 * @author Dahua Lin
 */


#include "distr_test_base.h"
#include <light_mat/random/poisson_distr.h>

default_rand_stream rstream;
const index_t N = 200000;

SIMPLE_CASE( test_poisson_naive )
{
	const double mu = 3.2;
	poisson_distr<uint32_t, naive_> distr(mu);

	ASSERT_EQ( distr.mean(), mu );
	ASSERT_EQ( distr.var(), mu );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, 8, ptol );
}


AUTO_TPACK( test_poissond )
{
	ADD_SIMPLE_CASE( test_poisson_naive )
}


