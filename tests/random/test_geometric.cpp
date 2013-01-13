/**
 * @file test_geometric.cpp
 *
 * @brief Unit testing of geometric distribution
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/geometric_distr.h>

default_rand_stream rstream;
const index_t N = 200000;

SIMPLE_CASE( test_geometric_naive )
{
	const double p = 0.4;
	geometric_distr<uint32_t, naive_> distr(p);

	ASSERT_EQ( distr.p(), p );
	ASSERT_APPROX( distr.mean(), (1.0 - p) / p, 1.0e-15 );
	ASSERT_APPROX( distr.var(), (1.0 - p) / (p * p), 1.0e-15 );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, 5, ptol );
}


AUTO_TPACK( test_geometric )
{
	ADD_SIMPLE_CASE( test_geometric_naive )
}


