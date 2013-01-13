/**
 * @file test_gammad.cpp
 *
 * @brief Unit testing of gamma distribution
 *
 * @author Dahua Lin
 */


#include "distr_test_base.h"
#include <light_mat/random/gamma_distr.h>

default_rand_stream rstream;
const index_t N = 200000;

T_CASE( test_std_gamma_g1_basic )
{
	T alpha = T(2.4);
	std_gamma_distr<T, basic_> distr(alpha);

	ASSERT_EQ( distr.alpha(), alpha );
	ASSERT_EQ( distr.beta(), T(1) );
	ASSERT_EQ( distr.mean(), alpha );
	ASSERT_EQ( distr.var(), alpha );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0 / alpha;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_std_gamma_l1_basic )
{
	T alpha = T(0.4);
	std_gamma_distr<T, basic_> distr(alpha);

	ASSERT_EQ( distr.alpha(), alpha );
	ASSERT_EQ( distr.beta(), T(1) );
	ASSERT_EQ( distr.mean(), alpha );
	ASSERT_EQ( distr.var(), alpha );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0 / alpha;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}


T_CASE( test_gamma_g1_basic )
{
	T alpha = T(2.4);
	T beta = T(1.6);
	gamma_distr<T, basic_> distr(alpha, beta);

	ASSERT_EQ( distr.alpha(), alpha );
	ASSERT_EQ( distr.beta(), beta );
	ASSERT_EQ( distr.mean(), alpha * beta );
	ASSERT_APPROX( distr.var(), alpha * beta * beta, 1.0e-15 );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0 / alpha;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}

T_CASE( test_gamma_l1_basic )
{
	T alpha = T(0.4);
	T beta = T(1.6);
	gamma_distr<T, basic_> distr(alpha, beta);

	ASSERT_EQ( distr.alpha(), alpha );
	ASSERT_EQ( distr.beta(), beta );
	ASSERT_EQ( distr.mean(), alpha * beta );
	ASSERT_APPROX( distr.var(), alpha * beta * beta, 1.0e-15 );

	double tol_mean = get_mean_tol(distr, N);
	double kappa = 6.0 / alpha;
	double tol_var = get_var_tol(distr, N, kappa);

	test_real_rng(distr, rstream, N, tol_mean, tol_var);
}


AUTO_TPACK( test_std_gamma_basic )
{
	ADD_T_CASE( test_std_gamma_g1_basic, double )
	ADD_T_CASE( test_std_gamma_g1_basic, float )
	ADD_T_CASE( test_std_gamma_l1_basic, double )
	ADD_T_CASE( test_std_gamma_l1_basic, float )
}

AUTO_TPACK( test_gamma_basic )
{
	ADD_T_CASE( test_gamma_g1_basic, double )
	ADD_T_CASE( test_gamma_g1_basic, float )
	ADD_T_CASE( test_gamma_l1_basic, double )
	ADD_T_CASE( test_gamma_l1_basic, float )
}



