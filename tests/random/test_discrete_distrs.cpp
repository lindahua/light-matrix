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
	dense_col<double> expect_p(K);
	for (index_t k = 0; k < K; ++k)
	{
		expect_p[k] = distr.p(k);
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


SIMPLE_CASE( std_uniform_int, simple )
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

SIMPLE_CASE( uniform_int, simple )
{
	const index_t a = 2;
	const index_t b = 6;
	std_uniform_int_distr<> distr(a, b);

	const index_t s = b - a + 1;

	ASSERT_EQ( distr.a(), 0 );
	ASSERT_EQ( distr.b(), b );
	ASSERT_EQ( distr.span(), s );
	ASSERT_EQ( distr.mean(), double(b) / 2 );
	ASSERT_EQ( distr.var(),  double((s * s - 1)) / 12 );

	test_discrete_rng(distr, rstream, N, b+1, default_ptol);
}


BEGIN_TPACK( uniform_int )
	ADD_SIMPLE_CASE( std_uniform_int, simple )
	ADD_SIMPLE_CASE( uniform_int, simple )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( uniform_int )
END_MAIN_SUITE




