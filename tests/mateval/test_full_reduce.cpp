/**
 * @file test_full_reduce.cpp
 *
 * Test of full reduction on matrices
 * 
 * @author Dahua Lin 
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_reduce.h>
#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

typedef atags::simd<default_simd_kind> simd_tag;

inline double randunif()
{
	double u = (double)std::rand() / double(RAND_MAX);
	return u * 2.0 - 1.0;
}

template<class Mat, typename T>
void fill_rand(IRegularMatrix<Mat, T>& mat)
{
	for (index_t j = 0; j < mat.ncolumns(); ++j)
	{
		for (index_t i = 0; i < mat.nrows(); ++i)
		{
			mat(i, j) = randunif();
		}
	}
}

const index_t max_len = 48;

#define DEF_FULL_REDUC_CASE(Name, UpdateExpr, initk, emptyval, tol ) \
		SIMPLE_CASE( full_reduce, Name ) { \
			dense_col<double> s(max_len); \
			fill_rand(s); \
			for (index_t k = initk; k <= max_len; ++k) { \
				auto sk = s(range(0, k)); \
				double r0 = emptyval; \
				for (index_t i = 0; i < k; ++i) UpdateExpr; \
				double r = Name(sk); \
				ASSERT_APPROX(r, r0, tol); \
				double r1 = Name(sk, linear_macc<atags::scalar>()); \
				ASSERT_APPROX(r1, r0, tol); \
				double r2 = Name(sk, linear_macc<simd_tag>()); \
				ASSERT_APPROX(r2, r0, tol); \
				double r3 = Name(sk, percol_macc<atags::scalar>()); \
				ASSERT_APPROX(r3, r0, tol); \
				double r4 = Name(sk, percol_macc<simd_tag>()); \
				ASSERT_APPROX(r4, r0, tol); } }

#define DEF_FULL_REDUC_CASE_2(Name, UpdateExpr, initk, emptyval, tol ) \
		SIMPLE_CASE( full_reduce, Name ) { \
			dense_col<double> s1(max_len); \
			dense_col<double> s2(max_len); \
			fill_rand(s1); \
			fill_rand(s2); \
			for (index_t k = initk; k <= max_len; ++k) { \
				auto sk1 = s1(range(0, k)); \
				auto sk2 = s2(range(0, k)); \
				double r0 = emptyval; \
				for (index_t i = 0; i < k; ++i) UpdateExpr; \
				double r = Name(sk1, sk2); \
				ASSERT_APPROX(r, r0, tol); \
				double r1 = Name(sk1, sk2, linear_macc<atags::scalar>()); \
				ASSERT_APPROX(r1, r0, tol); \
				double r2 = Name(sk1, sk2, linear_macc<simd_tag>()); \
				ASSERT_APPROX(r2, r0, tol); \
				double r3 = Name(sk1, sk2, percol_macc<atags::scalar>()); \
				ASSERT_APPROX(r3, r0, tol); \
				double r4 = Name(sk1, sk2, percol_macc<simd_tag>()); \
				ASSERT_APPROX(r4, r0, tol); } }



DEF_FULL_REDUC_CASE( sum, r0 += s[i], 0, 0.0, 1.0e-12 )
DEF_FULL_REDUC_CASE( maximum, r0 = math::max(r0, s[i]), 1, -1000.0, 1.0e-16 )
DEF_FULL_REDUC_CASE( minimum, r0 = math::min(r0, s[i]), 1,  1000.0, 1.0e-16 )

SIMPLE_CASE( full_reduce, mean )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		auto sk = s(range(0, k));

		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += s[i];
		r0 /= k;

		double r = mean(sk);
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = mean(sk, linear_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = mean(sk, linear_macc<simd_tag>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = mean(sk, percol_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = mean(sk, percol_macc<simd_tag>());
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


DEF_FULL_REDUC_CASE( asum, r0 += math::abs(s[i]), 0, 0.0, 1.0e-12 )
DEF_FULL_REDUC_CASE( amax, r0 = math::max(r0, math::abs(s[i])), 0, 0.0, 1.0e-16 )
DEF_FULL_REDUC_CASE( sqsum, r0 += math::sqr(s[i]), 0, 0.0, 1.0e-12 )

SIMPLE_CASE( full_reduce, amean )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		auto sk = s(range(0, k));

		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(s[i]);
		r0 /= k;

		double r = amean(sk);
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = amean(sk, linear_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = amean(sk, linear_macc<simd_tag>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = amean(sk, percol_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, 1.0e-12);

		r = amean(sk, percol_macc<simd_tag>());
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


DEF_FULL_REDUC_CASE_2( diff_asum, r0 += math::abs(sk1[i] - sk2[i]), 0, 0.0, 1.0e-12 )
DEF_FULL_REDUC_CASE_2( diff_amax, r0 = math::max(r0, math::abs(sk1[i] - sk2[i])), 0, 0.0, 1.0e-16 )
DEF_FULL_REDUC_CASE_2( diff_sqsum, r0 += math::sqr(sk1[i] - sk2[i]), 0, 0.0, 1.0e-12 )

SIMPLE_CASE( full_reduce, diff_amean )
{
	dense_col<double> s1(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s1);
	fill_rand(s2);

	double tol = 1.0e-12;

	for (index_t k = 1; k <= max_len; ++k)
	{
		auto sk1 = s1(range(0, k));
		auto sk2 = s2(range(0, k));

		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(sk1[i] - sk2[i]);
		r0 /= k;

		double r = diff_amean(sk1, sk2);
		ASSERT_APPROX(r, r0, tol);

		r = diff_amean(sk1, sk2, linear_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, tol);

		r = diff_amean(sk1, sk2, linear_macc<simd_tag>());
		ASSERT_APPROX(r, r0, tol);

		r = diff_amean(sk1, sk2, percol_macc<atags::scalar>());
		ASSERT_APPROX(r, r0, tol);

		r = diff_amean(sk1, sk2, percol_macc<simd_tag>());
		ASSERT_APPROX(r, r0, tol);
	}
}


DEF_FULL_REDUC_CASE_2( dot, r0 += sk1[i] * sk2[i], 0, 0.0, 1.0e-12 )


BEGIN_TPACK( full_reduce )
	ADD_SIMPLE_CASE( full_reduce, sum )
	ADD_SIMPLE_CASE( full_reduce, maximum )
	ADD_SIMPLE_CASE( full_reduce, minimum )
	ADD_SIMPLE_CASE( full_reduce, mean )

	ADD_SIMPLE_CASE( full_reduce, asum )
	ADD_SIMPLE_CASE( full_reduce, amax )
	ADD_SIMPLE_CASE( full_reduce, sqsum )
	ADD_SIMPLE_CASE( full_reduce, amean )

	ADD_SIMPLE_CASE( full_reduce, diff_asum )
	ADD_SIMPLE_CASE( full_reduce, diff_amax )
	ADD_SIMPLE_CASE( full_reduce, diff_sqsum )
	ADD_SIMPLE_CASE( full_reduce, diff_amean )

	ADD_SIMPLE_CASE( full_reduce, dot )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( full_reduce )
END_MAIN_SUITE

