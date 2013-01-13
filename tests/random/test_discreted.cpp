/**
 * @file test_discrete_distrs.cpp
 *
 * Unit testing of discrete distribution
 * 
 * @author Dahua Lin 
 */

#include "distr_test_base.h"
#include <light_mat/random/discrete_distr.h>

default_rand_stream rstream;
const index_t N = 200000;


SIMPLE_CASE( test_discrete_naive )
{
	discrete_distr<uint32_t, naive_> distr { 0.3, 0.6, 0.1, 0.7, 0.3 };

	ASSERT_EQ( distr.n(), 5 );
	ASSERT_APPROX( distr.p(0), 0.15, 1.0e-15 );
	ASSERT_APPROX( distr.p(1), 0.30, 1.0e-15 );
	ASSERT_APPROX( distr.p(2), 0.05, 1.0e-15 );
	ASSERT_APPROX( distr.p(3), 0.35, 1.0e-15 );
	ASSERT_APPROX( distr.p(4), 0.15, 1.0e-15 );
	ASSERT_EQ( distr.p(5), 0.00 );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, 6, ptol );
}

AUTO_TPACK( test_discreted )
{
	ADD_SIMPLE_CASE( test_discrete_naive )
}



