/**
 * @file test_bernoulli.cpp
 *
 * @brief Unit testing of bernoulli distribution
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/bernoulli_distr.h>

default_rand_stream rstream;
const index_t N = 200000;

SIMPLE_CASE( test_std_bernoulli )
{
	std_bernoulli_distr distr;

	ASSERT_EQ( distr.p(), 0.5 );
	ASSERT_EQ( distr.p(true), 0.5 );
	ASSERT_EQ( distr.p(false), 0.5 );
	ASSERT_EQ( distr.mean(), 0.5 );
	ASSERT_EQ( distr.var(), 0.25 );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, 2, ptol );
}

SIMPLE_CASE( test_bernoulli )
{
	const double p = 0.3;
	bernoulli_distr distr(p);

	ASSERT_EQ( distr.p(), p );
	ASSERT_EQ( distr.p(true), p );
	ASSERT_EQ( distr.p(false), 1.0 - p );
	ASSERT_EQ( distr.mean(), p );
	ASSERT_EQ( distr.var(), p * (1.0 - p) );

	double ptol = get_p_tol(N);
	test_discrete_rng(distr, rstream, N, 2, ptol );
}

AUTO_TPACK( test_bernoulli )
{
	ADD_SIMPLE_CASE( test_std_bernoulli )
	ADD_SIMPLE_CASE( test_bernoulli )
}
