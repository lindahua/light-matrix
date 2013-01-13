/**
 * @file test_real_distrs.cpp
 *
 * @brief Unit testing of real-valued distributions
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
void test_real_rng(const Distr& distr, RStream& rs, index_t n,
		double tol_mean, double tol_var, bool print_stats = false)
{
	static_assert(std::is_floating_point<typename Distr::result_type>::value,
			"Incorrect result_type.");

	double expect_mean = distr.mean();
	double expect_var = distr.var();

	double sx = 0.0;
	double sx2 = 0.0;

	for (index_t i = 0; i < n; ++i)
	{
		double x = distr(rs);
		sx += x;
		sx2 += math::sqr(x);
	}

	double actual_mean = sx / double(n);
	double actual_var = sx2 / double(n) - math::sqr(actual_mean);

	if (print_stats)
	{
		std::printf("\n");
		std::printf("amean = %g, emean = %g, tol = %g\n", actual_mean, expect_mean, tol_mean);
		std::printf("avar = %g, evar = %g, tol = %g\n", actual_var, expect_var, tol_var);
	}

	ASSERT_APPROX(actual_mean, expect_mean, tol_mean);
	ASSERT_APPROX(actual_var, expect_var, tol_var);
}


template<class Distr, class RStream, typename Kind>
void test_real_rng_simd(const Distr& distr, RStream& rs, Kind,
		index_t n, double tol_mean, double tol_var, bool print_stats=false)
{
	typedef typename Distr::result_type T;
	const unsigned int W = simd_traits<T, Kind>::pack_width;
	typedef simd_pack<T, Kind> pack_t;

	double expect_mean = distr.mean();
	double expect_var = distr.var();

	auto distr_ = simdize_map<Distr, Kind>::get(distr);

	double sx = 0.0;
	double sx2 = 0.0;

	LMAT_ALIGN(32) T sa[W];

	index_t m = n / (index_t)W;
	n = m * (index_t)W;

	for (index_t i = 0; i < m; ++i)
	{
		pack_t pk = distr_(rs);

		pk.store_a(sa);

		double cx = 0.0;
		double cx2 = 0.0;

		for (unsigned int j = 0; j < W; ++j)
		{
			double x = sa[j];
			cx += x;
			cx2 += math::sqr(x);
		}

		sx += cx;
		sx2 += cx2;
	}

	double actual_mean = sx / double(n);
	double actual_var = sx2 / double(n) - math::sqr(actual_mean);

	if (print_stats)
	{
		std::printf("\n");
		std::printf("amean = %g, emean = %g, tol = %g\n", actual_mean, expect_mean, tol_mean);
		std::printf("avar = %g, evar = %g, tol = %g\n", actual_var, expect_var, tol_var);
	}

	ASSERT_APPROX(actual_mean, expect_mean, tol_mean);
	ASSERT_APPROX(actual_var, expect_var, tol_var);
}


default_rand_stream rstream;
const index_t N = 200000;


template<class Distr>
inline double get_mean_tol(const Distr& distr, index_t n, double ratio=8.0)
{
	return ratio * math::sqrt(distr.var() / double(n));
}

template<class Distr>
inline double get_var_tol(const Distr& distr, index_t n, double kappa, double ratio=10.0)
{
	// kappa is excessive kurtosis

	return ratio * math::sqrt( (kappa + 2.0) * math::sqr(distr.var()) / double(n) );
}


/************************************************
 *
 *  Uniform real
 *
 ************************************************/

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


/************************************************
 *
 *  Exponential
 *
 ************************************************/

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


