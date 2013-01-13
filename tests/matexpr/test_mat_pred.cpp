/**
 * @file test_mat_pred.cpp
 *
 * @brief Unit testing of element-wise predicates on matrices
 *
 * @author Dahua Lin
 */

#include "matfun_test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_pred.h>
#include <light_mat/matexpr/mat_arith.h>
#include <limits>

#include <light_mat/matexpr/map_expr_inspect.h>

using namespace lmat;
using namespace lmat::test;

typedef mask_t<double> fmsk;


typedef std::numeric_limits<double> nlim;

const index_t n_spec_values = 8;
static double spec_values[n_spec_values] = {};

void fill_spec_values()
{
	spec_values[0] = nlim::infinity();
	spec_values[1] = - nlim::infinity();
	spec_values[2] = nlim::quiet_NaN();
	spec_values[3] = - nlim::quiet_NaN();
	spec_values[4] = 0.0;
	spec_values[5] = -0.0;
	spec_values[6] = 2.0;
	spec_values[7] = -2.0;
}


#define DEF_NUMCOMP_CASE( name, op, aexpr, bexpr, cexpr ) \
	MN_CASE( mat_##name ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		typedef dense_matrix<bool, M, N> bmat_t; \
		typedef dense_matrix<fmsk, M, N> mmat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); \
		mat_t B(m, n); \
		double c = double( cexpr ); \
		for (index_t i = 0; i < m * n; ++i) A[i] = double( aexpr ); \
		for (index_t i = 0; i < m * n; ++i) B[i] = double( bexpr ); \
		bmat_t AB_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) AB_r[i] = (A[i] op B[i]); \
		bmat_t AC_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) AC_r[i] = (A[i] op c); \
		bmat_t CB_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) CB_r[i] = (c op B[i]); \
		mmat_t ABm = (A op B); \
		bmat_t AB = to_bool(ABm); \
		ASSERT_EQ( AB.nrows(), m ); \
		ASSERT_EQ( AB.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, AB, AB_r ); \
		mmat_t ACm = (A op c); \
		bmat_t AC = to_bool(ACm); \
		ASSERT_EQ( AC.nrows(), m ); \
		ASSERT_EQ( AC.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, AC, AC_r ); \
		mmat_t CBm = (c op B); \
		bmat_t CB = to_bool(CBm); \
		ASSERT_EQ( CB.nrows(), m ); \
		ASSERT_EQ( CB.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, CB, CB_r ); }

#define DEF_NUMCOMP_TESTS( name, op ) \
	DEF_NUMCOMP_CASE( name, op, i+1, (i+1) * (i%2==1), (m*n)/2 ) \
	AUTO_TPACK( mat_##name ) { \
		ADD_MN_CASE_3X3( mat_##name, DM, DN ) \
	}


#define DEF_NUMPRED_TESTS( name ) \
	MN_CASE( mat_##name ) { \
		typedef dense_matrix<double, M, N> mat_t; \
		typedef dense_matrix<bool, M, N> bmat_t; \
		typedef dense_matrix<fmsk, M, N> mmat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		mat_t A(m, n); \
		fill_spec_values(); \
		for (index_t i = 0; i < m * n; ++i) A[i] = spec_values[i % n_spec_values]; \
		mmat_t Rm = name(A); \
		bmat_t R = to_bool(Rm); \
		bmat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) R_r[i] = math::name(A[i]); \
		ASSERT_EQ( R.nrows(), m ); \
		ASSERT_EQ( R.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R, R_r); } \
	AUTO_TPACK( mat_##name ) { \
		ADD_MN_CASE_3X3( mat_##name, DM, DN ) \
	}

#define DEF_LOGICAL_TESTS_1( name, op, bop ) \
	MN_CASE( mat_##name ) { \
		typedef dense_matrix<bool, M, N> bmat_t; \
		typedef dense_matrix<fmsk, M, N> mmat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		bmat_t A(m, n); \
		for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2); \
		mmat_t Am = to_f64m(A); \
		bmat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) R_r[i] = bop(A[i]); \
		bmat_t R = op(A); \
		mmat_t Rm = op(Am); \
		bmat_t Rmb = to_bool(Rm); \
		ASSERT_EQ( R.nrows(), m ); \
		ASSERT_EQ( R.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R, R_r); \
		ASSERT_EQ( Rmb.nrows(), m ); \
		ASSERT_EQ( Rmb.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, Rmb, R_r); } \
	AUTO_TPACK( mat_##name ) { \
		ADD_MN_CASE_3X3( mat_##name, DM, DN ) \
	}

#define DEF_LOGICAL_TESTS_2( name, op, bop ) \
	MN_CASE( mat_##name ) { \
		typedef dense_matrix<bool, M, N> bmat_t; \
		typedef dense_matrix<fmsk, M, N> mmat_t; \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		bmat_t A(m, n); \
		bmat_t B(m, n); \
		for (index_t i = 0; i < m * n; ++i) A[i] = bool(i % 2); \
		for (index_t i = 0; i < m * n; ++i) B[i] = bool(i % 3); \
		mmat_t Am = to_f64m(A); \
		mmat_t Bm = to_f64m(B); \
		bmat_t R_r(m, n); \
		for (index_t i = 0; i < m * n; ++i) R_r[i] = (A[i] bop B[i]); \
		bmat_t R1 = A op B; \
		mmat_t Rm2 = Am op Bm; \
		bmat_t R2 = to_bool(Rm2); \
		bmat_t R3 = Am op B; \
		bmat_t R4 = A op Bm; \
		ASSERT_EQ( R1.nrows(), m ); \
		ASSERT_EQ( R1.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R1, R_r); \
		ASSERT_EQ( R2.nrows(), m ); \
		ASSERT_EQ( R2.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R2, R_r); \
		ASSERT_EQ( R3.nrows(), m ); \
		ASSERT_EQ( R3.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R3, R_r); \
		ASSERT_EQ( R4.nrows(), m ); \
		ASSERT_EQ( R4.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, R1, R_r); \
    } \
	AUTO_TPACK( mat_##name ) { \
		ADD_MN_CASE_3X3( mat_##name, DM, DN ) \
	}


// specific tests

DEF_NUMCOMP_TESTS( eq, == )
DEF_NUMCOMP_TESTS( ne, != )
DEF_NUMCOMP_TESTS( ge, >= )
DEF_NUMCOMP_TESTS( gt, >  )
DEF_NUMCOMP_TESTS( le, <= )
DEF_NUMCOMP_TESTS( lt, <  )

DEF_NUMPRED_TESTS( signbit )
DEF_NUMPRED_TESTS( isfinite )
DEF_NUMPRED_TESTS( isinf )
DEF_NUMPRED_TESTS( isnan )

// ewise logical tests

DEF_LOGICAL_TESTS_1( logical_not, ~, !)
DEF_LOGICAL_TESTS_2( logical_and, &, && )
DEF_LOGICAL_TESTS_2( logical_or,  |, || )
DEF_LOGICAL_TESTS_2( logical_eq,  ==, == )
DEF_LOGICAL_TESTS_2( logical_ne,  !=, != )

// conditional selection

MN_CASE( mat_cond_b )
{
	typedef dense_matrix<bool, M, N> bmat_t;
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M; \
	const index_t n = N == 0 ? DN : N;

	bmat_t C(m, n);
	for (index_t i = 0; i < m * n; ++i) C[i] = bool(i % 2);

	mat_t X(m, n); fill_ran(X, 0.0, 1.0);
	mat_t Y(m, n); fill_ran(Y, 0.0, 1.0);
	double cv = X[0] + Y[0];

	typedef map_expr<ftags::cond_, bmat_t, mat_t, mat_t> expr_t;
	typedef preferred_macc_policy<expr_t> pmap;

	ASSERT_TRUE( pmap::prefer_linear );
	ASSERT_FALSE( pmap::prefer_simd );

	mat_t R1 = cond(C, X, Y);
	mat_t R1_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R1_r[i] = C[i] ? X[i] : Y[i];
	ASSERT_EQ( R1.nrows(), m );
	ASSERT_EQ( R1.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R1, R1_r );

	mat_t R2 = cond(C, X, cv);
	mat_t R2_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R2_r[i] = C[i] ? X[i] : cv;
	ASSERT_EQ( R2.nrows(), m );
	ASSERT_EQ( R2.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R2, R2_r );

	mat_t R3 = cond(C, cv, Y);
	mat_t R3_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R3_r[i] = C[i] ? cv : Y[i];
	ASSERT_EQ( R3.nrows(), m );
	ASSERT_EQ( R3.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R3, R3_r );
}

MN_CASE( mat_cond_m )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M; \
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = double(i + 1);
		B[i] = A[i] + double(i % 2);
	}

	mat_t X(m, n); fill_ran(X, 0.0, 1.0);
	mat_t Y(m, n); fill_ran(Y, 0.0, 1.0);
	double cv = X[0] + Y[0];

	typedef map_expr<ftags::cond_, map_expr<ftags::eq_, mat_t, mat_t>, mat_t, mat_t> expr_t;
	typedef preferred_macc_policy<expr_t> pmap;

	const int pw = simd_traits<double, default_simd_kind>::pack_width;
	ASSERT_TRUE( pmap::prefer_linear );
	bool use_simd = ((M * N) % pw == 0);

	ASSERT_EQ( pmap::prefer_simd, use_simd );

	mat_t R1 = cond(A == B, X, Y);
	mat_t R1_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R1_r[i] = A[i] == B[i] ? X[i] : Y[i];
	ASSERT_EQ( R1.nrows(), m );
	ASSERT_EQ( R1.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R1, R1_r );

	mat_t R2 = cond(A == B, X, cv);
	mat_t R2_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R2_r[i] = A[i] == B[i] ? X[i] : cv;
	ASSERT_EQ( R2.nrows(), m );
	ASSERT_EQ( R2.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R2, R2_r );

	mat_t R3 = cond(A == B, cv, Y);
	mat_t R3_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R3_r[i] = A[i] == B[i] ? cv : Y[i];
	ASSERT_EQ( R3.nrows(), m );
	ASSERT_EQ( R3.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, R3, R3_r );
}


AUTO_TPACK( mat_cond_b )
{
	ADD_MN_CASE_3X3( mat_cond_b, DM, DN )
}

AUTO_TPACK( mat_cond_m )
{
	ADD_MN_CASE_3X3( mat_cond_m, DM, DN )
}


