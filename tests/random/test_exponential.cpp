/**
 * @file test_exponential.cpp
 *
 * @brief Unit testing of exponential distributions
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/exponential_distr.h>

default_rand_stream rstream;
const index_t N = 200000;


T_CASE( test_std_exponential )
{
	std_exponential_distr<T> distr;

	ASSERT_EQ( distr.lambda(), T(1) );
	ASSERT_EQ( distr.beta(), T(1) );
	ASSERT_EQ( distr.mean(), T(1) );
	ASSERT_EQ( distr.var(), T(1) );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_std_exponential_sse )
{
	std_exponential_distr<T> distr;

	static_assert(is_simdizable<std_exponential_distr<T>, sse_t>::value,
			"std_exponential_distr should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_std_exponential_avx )
{
	std_exponential_distr<T> distr;

	static_assert(is_simdizable<std_exponential_distr<T>, avx_t>::value,
			"std_exponential_distr should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, avx_t(), N, tol_mean, tol_var);
}
#endif

T_CASE( test_exponential )
{
	T lambda = T(2);
	exponential_distr<T> distr(lambda);

	ASSERT_EQ( distr.lambda(), T(2) );
	ASSERT_EQ( distr.beta(), T(0.5) );
	ASSERT_EQ( distr.mean(), T(0.5) );
	ASSERT_EQ( distr.var(), T(0.25) );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}


T_CASE( test_exponential_sse )
{
	T lambda = T(2);
	exponential_distr<T> distr(lambda);

	static_assert(is_simdizable<exponential_distr<T>, sse_t>::value,
			"exponential_distr should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_exponential_avx )
{
	T lambda = T(2);
	exponential_distr<T> distr(lambda);

	static_assert(is_simdizable<exponential_distr<T>, avx_t>::value,
			"exponential_distr should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = -1.2;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, avx_t(), N, tol_mean, tol_var);
}
#endif


AUTO_TPACK( exponential )
{
	ADD_T_CASE( test_std_exponential, double )
	ADD_T_CASE( test_std_exponential, float )
	ADD_T_CASE( test_std_exponential_sse, double )
	ADD_T_CASE( test_std_exponential_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_std_exponential_avx, double )
	ADD_T_CASE( test_std_exponential_avx, float )
#endif

	ADD_T_CASE( test_exponential, double )
	ADD_T_CASE( test_exponential, float )
	ADD_T_CASE( test_exponential_sse, double )
	ADD_T_CASE( test_exponential_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_exponential_avx, double )
	ADD_T_CASE( test_exponential_avx, float )
#endif
}


