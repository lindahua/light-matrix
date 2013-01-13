/**
 * @file test_normal.cpp
 *
 * @brief Unit testing of normal distributions
 *
 * @author Dahua Lin
 */

#include "distr_test_base.h"
#include <light_mat/random/normal_distr.h>

default_rand_stream rstream(4321);
const index_t N = 200000;


T_CASE( test_std_normal_icdf )
{
	std_normal_distr<T, icdf_> distr;

	ASSERT_EQ( distr.mean(), T(0) );
	ASSERT_EQ( distr.stddev(), T(1) );
	ASSERT_EQ( distr.var(), T(1) );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_std_normal_icdf_sse )
{
	std_normal_distr<T, icdf_> distr;

	static_assert(is_simdizable<std_normal_distr<T, icdf_>, sse_t>::value,
			"std_normal_distr (icdf_) should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_std_normal_icdf_avx )
{
	std_normal_distr<T, icdf_> distr;

	static_assert(is_simdizable<std_normal_distr<T, icdf_>, avx_t>::value,
			"std_normal_distr (icdf_) should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}
#endif


T_CASE( test_normal_icdf_sse )
{
	normal_distr<T, icdf_> distr(T(2.5), T(2.0));

	static_assert(is_simdizable<normal_distr<T, icdf_>, sse_t>::value,
			"normal_distr (icdf_) should be simdizable with sse");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}

#ifdef LMAT_HAS_AVX
T_CASE( test_normal_icdf_avx )
{
	normal_distr<T, icdf_> distr(T(2.5), T(2.0));

	static_assert(is_simdizable<normal_distr<T, icdf_>, avx_t>::value,
			"normal_distr (icdf_) should be simdizable with avx");

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng_simd(distr, rstream, sse_t(), N, tol_mean, tol_var);
}
#endif


T_CASE( test_normal_icdf )
{
	T mu = T(2.5);
	T sigma = T(2.0);
	normal_distr<T> distr(mu, sigma);

	ASSERT_EQ( distr.mean(), mu );
	ASSERT_EQ( distr.stddev(), sigma );
	ASSERT_EQ( distr.var(), sigma * sigma );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 0.0;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}


AUTO_TPACK( test_normal_icdf )
{
	ADD_T_CASE( test_std_normal_icdf, double )
	ADD_T_CASE( test_std_normal_icdf, float )
	ADD_T_CASE( test_std_normal_icdf_sse, double )
	ADD_T_CASE( test_std_normal_icdf_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_std_normal_icdf_avx, double )
	ADD_T_CASE( test_std_normal_icdf_avx, float )
#endif

	ADD_T_CASE( test_normal_icdf, double )
	ADD_T_CASE( test_normal_icdf, float )
	ADD_T_CASE( test_normal_icdf_sse, double )
	ADD_T_CASE( test_normal_icdf_sse, float )
#ifdef LMAT_HAS_AVX
	ADD_T_CASE( test_normal_icdf_avx, double )
	ADD_T_CASE( test_normal_icdf_avx, float )
#endif
}

