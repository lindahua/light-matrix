/**
 * @file test_uniform_real.cpp
 *
 * @brief Test uniform real distributions
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/uniform_real_distr.h>

default_rand_stream rstream;
const index_t N = 200000;


T_CASE( test_std_uniform_real )
{
	std_uniform_real_distr<T> distr;

	ASSERT_EQ( distr.a(), T(0) );
	ASSERT_EQ( distr.b(), T(1) );
	ASSERT_EQ( distr.span(), T(1) );
	ASSERT_EQ( distr.mean(), T(0.5) );
	ASSERT_EQ( distr.var(), T(1.0/12) );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_std_uniform_real_sse )
{
	std_uniform_real_distr<T> distr;

	static_assert(is_simdizable<std_uniform_real_distr<T>, sse_t>::value,
			"std_uniform_real_distr should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_std_uniform_real_avx )
{
	std_uniform_real_distr<T> distr;

	static_assert(is_simdizable<std_uniform_real_distr<T>, avx_t>::value,
			"std_uniform_real_distr should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, avx_t(), N, tol_mean, tol_var);
}
#endif


T_CASE( test_uniform_real )
{
	T a = T(1.6);
	T b = T(3.2);
	uniform_real_distr<T> distr(a, b);

	ASSERT_EQ( distr.a(), a );
	ASSERT_EQ( distr.b(), b );
	ASSERT_EQ( distr.span(), b-a );
	ASSERT_EQ( distr.mean(), (a+b)/T(2) );
	ASSERT_APPROX( distr.var(), math::sqr(b-a)/T(12), 1.0e-6 );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_uniform_real_sse )
{
	uniform_real_distr<T> distr(T(1.6), T(3.2));

	static_assert(is_simdizable<uniform_real_distr<T>, sse_t>::value,
			"uniform_real_distr should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_uniform_real_avx )
{
	uniform_real_distr<T> distr(T(1.6), T(3.2));

	static_assert(is_simdizable<uniform_real_distr<T>, avx_t>::value,
			"uniform_real_distr should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, avx_t(), N, tol_mean, tol_var);
}
#endif


AUTO_TPACK( uniform_real )
{
	ADD_T_CASE( test_std_uniform_real, double )
	ADD_T_CASE( test_std_uniform_real, float )
	ADD_T_CASE( test_std_uniform_real_sse, double )
	ADD_T_CASE( test_std_uniform_real_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_std_uniform_real_avx, double )
	ADD_T_CASE( test_std_uniform_real_avx, float )
#endif

	ADD_T_CASE( test_uniform_real, double )
	ADD_T_CASE( test_uniform_real, float )
	ADD_T_CASE( test_uniform_real_sse, double )
	ADD_T_CASE( test_uniform_real_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_uniform_real_avx, double )
	ADD_T_CASE( test_uniform_real_avx, float )
#endif
}



