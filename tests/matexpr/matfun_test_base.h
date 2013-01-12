/**
 * @file matfun_test_base.h
 *
 * @brief Basic support of unit testing for matrix-map functions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATFUN_TEST_BASE_H_
#define LIGHTMAT_MATFUN_TEST_BASE_H_

#include "../test_base.h"

using namespace lmat;
using namespace lmat::test;

const int DM = 8;
const int DN = 6;
const index_t LDim = 12;

template<class Expr>
inline void check_policy(const Expr& expr)
{
	typedef preferred_macc_policy<Expr> pmap;
	typedef decltype(expr.tag()) tag_t;
	typedef default_simd_kind skind;

	typedef typename matrix_traits<Expr>::value_type T;
	const int pw = (int)simd_traits<T, skind>::pack_width;
	ASSERT_TRUE( pmap::prefer_linear );

	int M = meta::nrows<Expr>::value;
	int N = meta::ncols<Expr>::value;

	bool use_simd = ((M * N) % pw == 0) &&
			meta::has_simd_support<tag_t, T, skind>::value;
	ASSERT_EQ( pmap::prefer_simd, use_simd );
}

#define DEFINE_MATFUN_TESTS1( matfun, scafun, al, ah, tol ) \
	MN_CASE( mat_math, matfun ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); fill_ran(A, al, ah); \
		check_policy( matfun(A) ); \
		mat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) \
			R_r[i] = scafun(A[i]); \
		mat_t R = matfun(A); \
		ASSERT_TRUE( my_is_approx(R, R_r, tol) ); \
	} \
	BEGIN_TPACK( mat##_matfun ) \
		ADD_MN_CASE_3X3( mat_math, matfun, DM, DN ) \
	END_TPACK

#define DEFINE_MATFUN_TESTS2( matfun, scafun, al, ah, bl, bh, tol ) \
	MN_CASE( mat_math, scafun ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); fill_ran(A, al, ah); \
		check_policy( matfun(A) ); \
		mat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) \
			R_r[i] = scafun(A[i]); \
		mat_t R = matfun(A); \
		ASSERT_TRUE( my_is_approx(R, R_r, tol) ); \
	} \
	BEGIN_TPACK( mat##_matfun ) \
		ADD_MN_CASE_3X3( mat_math, matfun, DM, DN ) \
	END_TPACK

#endif
