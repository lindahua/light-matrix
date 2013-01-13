/**
 * @file test_discrete_distrs.cpp
 *
 * Test PRNGs of discrete random distributions
 * 
 * @author Dahua Lin 
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/random/sfmt.h>

#include <light_mat/random/distributions.h>

using namespace lmat;
using namespace lmat::random;
using namespace lmat::test;


template<class Distr, class RStream>
void test_discrete_rng(const Distr& distr, RStream& rs, index_t n, index_t K, double ptol)
{
	typedef typename Distr::result_type RT;

	dense_col<double> expect_p(K);
	for (index_t k = 0; k < K; ++k)
	{
		expect_p[k] = distr.p((RT)k);
	}

	dense_col<uint32_t> counts(K, zero());

	for (index_t i = 0; i < n; ++i)
	{
		index_t x = (index_t)(distr(rs));
		if (x >= 0 && x < K) ++counts[x];
	}

	dense_col<double> actual_p(K);
	for (index_t k = 0; k < K; ++k)
	{
		actual_p[k] = double(counts[k]) / double(n);
	}

	ASSERT_VEC_APPROX(K, actual_p, expect_p, ptol);
}


default_rand_stream rstream;
const index_t N = 200000;
const double default_ptol = 5.0 / std::sqrt(double(N));


SIMPLE_CASE( test_std_uniform_int )
{
	const index_t b = 5;
	std_uniform_int_distr<> distr(b);

	const index_t s = b + 1;

	ASSERT_EQ( distr.a(), 0 );
	ASSERT_EQ( distr.b(), b );
	ASSERT_EQ( distr.span(), s );
	ASSERT_EQ( distr.mean(), double(b) / 2 );
	ASSERT_EQ( distr.var(),  double((s * s - 1)) / 12 );

	test_discrete_rng(distr, rstream, N, b+1, default_ptol);
}

SIMPLE_CASE( test_uniform_int )
{
	const index_t a = 2;
	const index_t b = 6;
	uniform_int_distr<> distr(a, b);

	const index_t s = b - a + 1;

	ASSERT_EQ( distr.a(), a );
	ASSERT_EQ( distr.b(), b );
	ASSERT_EQ( distr.span(), s );
	ASSERT_EQ( distr.mean(), double(a + b) / 2 );
	ASSERT_EQ( distr.var(),  double((s * s - 1)) / 12 );

	test_discrete_rng(distr, rstream, N, b+1, default_ptol);
}


SIMPLE_CASE( test_std_bernoulli )
{
	std_bernoulli_distr distr;

	ASSERT_EQ( distr.p(), 0.5 );
	ASSERT_EQ( distr.p(true), 0.5 );
	ASSERT_EQ( distr.p(false), 0.5 );
	ASSERT_EQ( distr.mean(), 0.5 );
	ASSERT_EQ( distr.var(), 0.25 );

	test_discrete_rng(distr, rstream, N, 2, default_ptol );
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

	test_discrete_rng(distr, rstream, N, 2, default_ptol );
}

SIMPLE_CASE( test_binomial_naive )
{
	uint32_t t = 5;
	const double p = 0.4;
	binomial_distr<uint32_t, naive_> distr(t, p);

	ASSERT_EQ( distr.t(), t );
	ASSERT_EQ( distr.p(), p );
	ASSERT_EQ( distr.mean(), double(t) * p );
	ASSERT_APPROX( distr.var(), double(t) * p * (1-p), 1.0e-15 );

	test_discrete_rng(distr, rstream, N, (index_t)t+1, default_ptol );
}


SIMPLE_CASE( test_geometric_naive )
{
	const double p = 0.4;
	geometric_distr<uint32_t, naive_> distr(p);

	ASSERT_EQ( distr.p(), p );
	ASSERT_APPROX( distr.mean(), (1.0 - p) / p, 1.0e-15 );
	ASSERT_APPROX( distr.var(), (1.0 - p) / (p * p), 1.0e-15 );

	test_discrete_rng(distr, rstream, N, 5, default_ptol );
}


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

	test_discrete_rng(distr, rstream, N, 6, default_ptol );
}


AUTO_TPACK( uniform_int )
{
	ADD_SIMPLE_CASE( test_std_uniform_int )
	ADD_SIMPLE_CASE( test_uniform_int )
}

AUTO_TPACK( bernoulli )
{
	ADD_SIMPLE_CASE( test_std_bernoulli )
	ADD_SIMPLE_CASE( test_bernoulli )
}

AUTO_TPACK( binomial )
{
	ADD_SIMPLE_CASE( test_binomial_naive )
}

AUTO_TPACK( geometric )
{
	ADD_SIMPLE_CASE( test_geometric_naive )
}

AUTO_TPACK( discrete )
{
	ADD_SIMPLE_CASE( test_discrete_naive )
}



