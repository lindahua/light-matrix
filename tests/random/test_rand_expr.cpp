/**
 * @file test_rand_expr.cpp
 *
 * @brief Unit testing of random matrix expressions
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/random/rand_expr.h>

using namespace lmat;
using namespace lmat::random;
using namespace lmat::test;

default_rand_stream rstream;

const unsigned int seed = 4321;

const index_t DM = 6;
const index_t DN = 8;

template<class Distr, class RStream, int M, int N>
inline void check_policy(const rand_expr<Distr, RStream, M, N>& expr)
{
	typedef rand_expr<Distr, RStream, M, N> expr_t;
	bool expect_usimd = is_simdizable<Distr, default_simd_kind>::value;
	// bool expect_usimd = true;

	typedef preferred_macc_policy<expr_t> pmap;
	ASSERT_EQ( pmap::prefer_simd, expect_usimd );
}

template<class Distr, class RStream, int M, int N>
void test_rand_expr(const rand_expr<Distr, RStream, M, N>& expr, Distr& distr0)
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename Distr::result_type T;

	check_policy(expr);

	rstream.set_seed(seed);

	dense_matrix<T, M, N> R_r(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		R_r[i] = distr0(rstream);
	}

	rstream.set_seed(seed);

	dense_matrix<T, M, N> R = expr;

	T tol = sizeof(T) == 4 ? T(1.0e-6) : T(1.0e-12);

	ASSERT_EQ( R.nrows(), m );
	ASSERT_EQ( R.ncolumns(), n );
	ASSERT_MAT_APPROX( m, n, R, R_r, tol );
}


/************************************************
 *
 *  test cases for standard rand_mat functions
 *
 ************************************************/

TMN_CASE( test_rand_mat )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	std_uniform_real_distr<T> distr0;
	matrix_shape<M, N> shape(m, n);
	test_rand_expr(rand_mat(distr0, rstream, shape), distr0);
}

AUTO_TPACK( test_rand_mat )
{
	ADD_TMN_CASE_3X3( test_rand_mat, double, DM, DN )
	ADD_TMN_CASE_3X3( test_rand_mat, float, DM, DN )
}


/************************************************
 *
 *  test cases for convenient functions
 *
 ************************************************/

// randu

SIMPLE_CASE( test_randu )
{
	std_uniform_real_distr<double> distr0;
	test_rand_expr(randu(rstream, DM, DN), distr0);
}

SIMPLE_CASE( test_randuf )
{
	std_uniform_real_distr<float> distr0;
	test_rand_expr(randuf(rstream, DM, DN), distr0);
}

T_CASE( test_randux )
{
	T a = T(2);
	T b = T(5);
	uniform_real_distr<T> distr0(a, b);
	test_rand_expr(randu(rstream, DM, DN, a, b), distr0);
}

AUTO_TPACK( test_randu )
{
	ADD_SIMPLE_CASE( test_randu )
	ADD_SIMPLE_CASE( test_randu )
	ADD_T_CASE_FP( test_randux )
}

// randn

SIMPLE_CASE( test_randn )
{
	std_normal_distr<double> distr0;
	test_rand_expr(randn(rstream, DM, DN), distr0);
}

SIMPLE_CASE( test_randnf )
{
	std_normal_distr<float> distr0;
	test_rand_expr(randnf(rstream, DM, DN), distr0);
}

T_CASE( test_randnx )
{
	T mu = T(2);
	T sigma = T(1.5);
	normal_distr<T> distr0(mu, sigma);
	test_rand_expr(randn(rstream, DM, DN, mu, sigma), distr0);
}

AUTO_TPACK( test_randn )
{
	ADD_SIMPLE_CASE( test_randn )
	ADD_SIMPLE_CASE( test_randn )
	ADD_T_CASE_FP( test_randnx )
}


// rande

SIMPLE_CASE( test_rande )
{
	std_exponential_distr<double> distr0;
	test_rand_expr(rande(rstream, DM, DN), distr0);
}

SIMPLE_CASE( test_randef )
{
	std_exponential_distr<float> distr0;
	test_rand_expr(randef(rstream, DM, DN), distr0);
}

T_CASE( test_randex )
{
	T lam = T(2);
	exponential_distr<T> distr0(lam);
	test_rand_expr(rande(rstream, DM, DN, lam), distr0);
}

AUTO_TPACK( test_rande )
{
	ADD_SIMPLE_CASE( test_rande )
	ADD_SIMPLE_CASE( test_rande )
	ADD_T_CASE_FP( test_randex )
}



// randg


T_CASE( test_randgx )
{
	T a = T(2);
	std_gamma_distr<T> distr0(a);
	test_rand_expr(randg(rstream, DM, DN, a), distr0);
}


T_CASE( test_randgx2 )
{
	T a = T(2);
	T b = T(1.5);
	gamma_distr<T> distr0(a, b);
	test_rand_expr(randg(rstream, DM, DN, a, b), distr0);
}


AUTO_TPACK( test_randg )
{
	ADD_T_CASE_FP( test_randgx )
	ADD_T_CASE_FP( test_randgx2 )
}




