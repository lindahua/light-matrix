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

#include <light_mat/mateval/ewise_eval.h>

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


template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}


#define DEFINE_MATFUN_TESTS1( fun, al, ah, tol ) \
	MN_CASE( mat_##fun ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); fill_ran(A, al, ah); \
		check_policy( fun(A) ); \
		mat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) \
			R_r[i] = math::fun(A[i]); \
		mat_t R = fun(A); \
		ASSERT_EQ( R.nrows(), m ); \
		ASSERT_EQ( R.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, R, R_r, tol ); \
	} \
	AUTO_TPACK( mat_##fun ) { \
		ADD_MN_CASE_3X3( mat_##fun, DM, DN ) \
	}

#define DEFINE_MATFUN_TESTS2( fun, al, ah, bl, bh, tol ) \
	MN_CASE( mat_##fun ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); fill_ran(A, al, ah); \
		mat_t B(m ,n); fill_ran(B, bl, bh); \
		double ca = (al + ah) * 0.5; \
		double cb = (bl + bh) * 0.5; \
		check_policy( fun(A, B) ); \
		check_policy( fun(A, cb) ); \
		check_policy( fun(ca, B) ); \
		mat_t AB_r(m, n); \
		mat_t AC_r(m, n); \
		mat_t CB_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) { \
			AB_r[i] = math::fun(A[i], B[i]); \
			AC_r[i] = math::fun(A[i], cb); \
			CB_r[i] = math::fun(ca, B[i]); \
		} \
		mat_t AB = fun(A, B); \
		mat_t AC = fun(A, cb); \
		mat_t CB = fun(ca, B); \
		ASSERT_EQ( AB.nrows(), m ); \
		ASSERT_EQ( AB.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, AB, AB_r, tol ); \
		ASSERT_EQ( AC.nrows(), m ); \
		ASSERT_EQ( AC.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, AC, AC_r, tol ); \
		ASSERT_EQ( CB.nrows(), m ); \
		ASSERT_EQ( CB.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, CB, CB_r, tol ); \
	} \
	AUTO_TPACK( mat_##fun ) { \
		ADD_MN_CASE_3X3( mat_##fun, DM, DN ) \
	}


#define DEFINE_MATFUN_TESTS3a( fun, vl, vh, tol ) \
	MN_CASE( mat_##fun ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); fill_ran(A, vl, vh); \
		mat_t B(m ,n); fill_ran(B, vl, vh); \
		mat_t C(m, n); fill_ran(C, vl, vh); \
		double a = A[0]; \
		double b = B[0]; \
		double c = C[0]; \
		check_policy( fun(A, B, C) ); \
		check_policy( fun(a, B, C) ); \
		check_policy( fun(A, b, C) ); \
		check_policy( fun(A, B, c) ); \
		check_policy( fun(a, b, C) ); \
		check_policy( fun(a, B, c) ); \
		check_policy( fun(A, b, c) ); \
		mat_t ABC_r(m, n); \
		mat_t aBC_r(m, n); \
		mat_t AbC_r(m, n); \
		mat_t ABc_r(m, n); \
		mat_t abC_r(m, n); \
		mat_t aBc_r(m, n); \
		mat_t Abc_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) { \
			ABC_r[i] = math::fun(A[i], B[i], C[i]); \
			aBC_r[i] = math::fun(a, B[i], C[i]); \
			AbC_r[i] = math::fun(A[i], b, C[i]); \
			ABc_r[i] = math::fun(A[i], B[i], c); \
			abC_r[i] = math::fun(a, b, C[i]); \
			aBc_r[i] = math::fun(a, B[i], c); \
			Abc_r[i] = math::fun(A[i], b, c); \
		} \
		mat_t ABC = fun(A, B, C); \
		mat_t aBC = fun(a, B, C); \
		mat_t AbC = fun(A, b, C); \
		mat_t ABc = fun(A, B, c); \
		mat_t abC = fun(a, b, C); \
		mat_t aBc = fun(a, B, c); \
		mat_t Abc = fun(A, b, c); \
		ASSERT_EQ( ABC.nrows(), m ); \
		ASSERT_EQ( ABC.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, ABC, ABC_r, tol ); \
		ASSERT_EQ( aBC.nrows(), m ); \
		ASSERT_EQ( aBC.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, aBC, aBC_r, tol ); \
		ASSERT_EQ( AbC.nrows(), m ); \
		ASSERT_EQ( AbC.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, AbC, AbC_r, tol ); \
		ASSERT_EQ( ABc.nrows(), m ); \
		ASSERT_EQ( ABc.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, ABc, ABc_r, tol ); \
		ASSERT_EQ( abC.nrows(), m ); \
		ASSERT_EQ( abC.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, abC, abC_r, tol ); \
		ASSERT_EQ( aBc.nrows(), m ); \
		ASSERT_EQ( aBc.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, aBc, aBc_r, tol ); \
		ASSERT_EQ( Abc.nrows(), m ); \
		ASSERT_EQ( Abc.ncolumns(), n); \
		ASSERT_MAT_APPROX( m, n, Abc, Abc_r, tol ); \
	} \
	AUTO_TPACK( mat_##fun ) { \
		ADD_MN_CASE_3X3( mat_##fun, DM, DN ) \
	}


#endif
