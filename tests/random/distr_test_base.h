/**
 * @file distr_test_base.h
 *
 * @brief Common facilities to help testing distributions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef DISTR_TEST_BASE_H_
#define DISTR_TEST_BASE_H_

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/math/math_base.h>
#include <light_mat/random/distr_fwd.h>

using namespace lmat;
using namespace lmat::test;
using namespace lmat::random;

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


inline double get_p_tol(index_t n)
{
	return 5.0 / std::sqrt(double(n));
}

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


#endif
