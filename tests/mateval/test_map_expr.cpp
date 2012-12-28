/**
 * @file test_map_expr.cpp
 *
 * @brief Unit testing of the framework for map expression
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#define DEFAULT_M_VALUE 9
#define DEFAULT_N_VALUE 8

#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/map_expr.h>

using namespace lmat;
using namespace lmat::test;

typedef atags::scalar scalar_tag;
typedef atags::simd<default_simd_kind> simd_tag;


template<class FTag, class A>
bool my_use_linear(const map_expr<FTag, A>& expr)
{
	const int M = meta::nrows<A>::value;
	const int N = meta::ncols<A>::value;

	const bool is_cont = meta::is_continuous<A>::value;

	return is_cont || M == 1 || N == 1;
}


template<class FTag, class A, class B>
bool my_use_linear(const map_expr<FTag, A, B>& expr)
{
	const int M = meta::common_nrows<A, B>::value;
	const int N = meta::common_ncols<A, B>::value;

	const bool is_cont =
			meta::is_continuous<A>::value &&
			meta::is_continuous<B>::value;

	return is_cont || M == 1 || N == 1;
}


template<class FTag, class A>
bool my_use_simd(const map_expr<FTag, A>& expr)
{
	bool use_linear = my_use_linear(expr);

	int pw = math::simd_traits<double, default_simd_kind>::pack_width;

	if (use_linear)
	{
		const int L = meta::nelems<A>::value;
		const bool is_cont = meta::is_continuous<A>::value;

		return is_cont && (L % pw == 0);
	}
	else
	{
		const int M = meta::nrows<A>::value;
		const bool is_cont_pc = meta::is_percol_continuous<A>::value;

		return is_cont_pc && (M % pw == 0);
	}
}

template<class FTag, class A, class B>
bool my_use_simd(const map_expr<FTag, A, B>& expr)
{
	bool use_linear = my_use_linear(expr);

	int pw = math::simd_traits<double, default_simd_kind>::pack_width;

	if (use_linear)
	{
		const int L = meta::common_nelems<A, B>::value;
		const bool is_cont =
				meta::is_continuous<A>::value &&
				meta::is_continuous<B>::value;

		return is_cont && (L % pw == 0);
	}
	else
	{
		const int M = meta::common_nrows<A, B>::value;
		const bool is_cont_pc =
				meta::is_percol_continuous<A>::value &&
				meta::is_percol_continuous<B>::value;

		return is_cont_pc && (M % pw == 0);
	}
}




template<typename STag1, typename DTag, int M, int N>
void test_mapexpr_1()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef mat_host<STag1, double, M, N> shost1_t;
	typedef mat_host<DTag, double, M, N> dhost_t;

	typedef typename shost1_t::cmat_t smat1_t;
	typedef typename dhost_t::mat_t dmat_t;

	shost1_t s1_h(m, n);
	dhost_t d_h(m, n);

	s1_h.fill_rand();

	smat1_t s1 = s1_h.get_cmat();
	dmat_t d = d_h.get_mat();

	typedef map_expr<sqr_, smat1_t> expr_t;
	expr_t e = make_map_expr(sqr_(), s1);

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n);

	// policy

	typedef preferred_macc_policy<expr_t> pmap;
	ASSERT_EQ( pmap::prefer_linear, my_use_linear(e) );
	ASSERT_EQ( pmap::prefer_simd, my_use_simd(e) );

	// evaluation

	d = e;

	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::sqr(s1(i, j));
	}

	double tol = 1.0e-14;
	ASSERT_MAT_APPROX(m, n, d, r, tol);
}


template<typename STag1, typename STag2, typename DTag, int M, int N>
void test_mapexpr_2()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef mat_host<STag1, double, M, N> shost1_t;
	typedef mat_host<STag2, double, M, N> shost2_t;
	typedef mat_host<DTag, double, M, N> dhost_t;

	typedef typename shost1_t::cmat_t smat1_t;
	typedef typename shost2_t::cmat_t smat2_t;
	typedef typename dhost_t::mat_t dmat_t;

	shost1_t s1_h(m, n);
	shost2_t s2_h(m, n);
	dhost_t d_h(m, n);

	s1_h.fill_rand();
	s2_h.fill_rand();

	smat1_t s1 = s1_h.get_cmat();
	smat2_t s2 = s2_h.get_cmat();
	dmat_t d = d_h.get_mat();

	typedef map_expr<sub_, smat1_t, smat2_t> expr_t;
	expr_t e = make_map_expr(sub_(), s1, s2);

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n);

	// policy

	typedef preferred_macc_policy<expr_t> pmap;
	ASSERT_EQ( pmap::prefer_linear, my_use_linear(e) );
	ASSERT_EQ( pmap::prefer_simd, my_use_simd(e) );

	// evaluation

	d = e;

	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = s1(i, j) - s2(i, j);
	}

	double tol = 1.0e-14;
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	// other forms

	double cv = 2.5;

	d = make_map_expr_fix2(sub_(), s1, cv);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = s1(i, j) - cv;
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix1(sub_(), cv, s2);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = cv - s2(i, j);
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

}


template<typename STag1, typename STag2, typename STag3, typename DTag, int M, int N>
void test_mapexpr_3()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef mat_host<STag1, double, M, N> shost1_t;
	typedef mat_host<STag2, double, M, N> shost2_t;
	typedef mat_host<STag3, double, M, N> shost3_t;
	typedef mat_host<DTag, double, M, N> dhost_t;

	typedef typename shost1_t::cmat_t smat1_t;
	typedef typename shost2_t::cmat_t smat2_t;
	typedef typename shost3_t::cmat_t smat3_t;
	typedef typename dhost_t::mat_t dmat_t;

	shost1_t s1_h(m, n);
	shost2_t s2_h(m, n);
	shost3_t s3_h(m, n);
	dhost_t d_h(m, n);

	s1_h.fill_rand();
	s2_h.fill_rand();
	s3_h.fill_rand();

	smat1_t s1 = s1_h.get_cmat();
	smat2_t s2 = s2_h.get_cmat();
	smat3_t s3 = s3_h.get_cmat();
	dmat_t d = d_h.get_mat();

	typedef map_expr<fma_, smat1_t, smat2_t, smat3_t> expr_t;
	expr_t e = make_map_expr(fma_(), s1, s2, s3);

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n);

	d = e;

	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = s1(i, j) * s2(i, j) + s3(i, j);
	}

	double tol = 1.0e-14;
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	// other forms

	double v1 = 2.5;
	double v2 = 3.2;
	double v3 = -1.8;

	d = make_map_expr_fix1(fma_(), v1, s2, s3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(v1, s2(i, j), s3(i, j));
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix2(fma_(), s1, v2, s3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(s1(i, j), v2, s3(i, j));
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix3(fma_(), s1, s2, v3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(s1(i, j), s2(i, j), v3);
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix12(fma_(), v1, v2, s3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(v1, v2, s3(i, j));
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix13(fma_(), v1, s2, v3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(v1, s2(i, j), v3);
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);

	d = make_map_expr_fix23(fma_(), s1, v2, v3);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = math::fma(s1(i, j), v2, v3);
	}
	ASSERT_MAT_APPROX(m, n, d, r, tol);
}



// Unary expressions

#define DEF_MEXPR_TESTS_1( stag, dtag ) \
		MN_CASE( map_expr, unary_##stag##_##dtag ) { test_mapexpr_1<stag, dtag, M, N>(); } \
		BEGIN_TPACK( unary_map_expr_##stag##_##dtag ) \
		ADD_MN_CASE_3X3( map_expr, unary_##stag##_##dtag, DM, DN ) \
		END_TPACK

DEF_MEXPR_TESTS_1( cont, cont )
DEF_MEXPR_TESTS_1( cont, bloc )
DEF_MEXPR_TESTS_1( cont, grid )
DEF_MEXPR_TESTS_1( bloc, cont )
DEF_MEXPR_TESTS_1( bloc, bloc )
DEF_MEXPR_TESTS_1( bloc, grid )
DEF_MEXPR_TESTS_1( grid, cont )
DEF_MEXPR_TESTS_1( grid, bloc )
DEF_MEXPR_TESTS_1( grid, grid )

// Binary expression

#define DEF_MEXPR_TESTS_2( stag1, stag2, dtag ) \
		MN_CASE( map_expr, binary_##stag1##_##stag2##_##dtag ) { test_mapexpr_2<stag1, stag2, dtag, M, N>(); } \
		BEGIN_TPACK( binary_map_expr_##stag1##_##stag2##_##dtag ) \
		ADD_MN_CASE_3X3( map_expr, binary_##stag1##_##stag2##_##dtag, DM, DN ) \
		END_TPACK

DEF_MEXPR_TESTS_2( cont, cont, cont )
DEF_MEXPR_TESTS_2( cont, cont, bloc )
DEF_MEXPR_TESTS_2( cont, cont, grid )
DEF_MEXPR_TESTS_2( cont, bloc, cont )
DEF_MEXPR_TESTS_2( cont, bloc, bloc )
DEF_MEXPR_TESTS_2( cont, bloc, grid )
DEF_MEXPR_TESTS_2( cont, grid, cont )
DEF_MEXPR_TESTS_2( cont, grid, bloc )
DEF_MEXPR_TESTS_2( cont, grid, grid )

DEF_MEXPR_TESTS_2( bloc, cont, cont )
DEF_MEXPR_TESTS_2( bloc, cont, bloc )
DEF_MEXPR_TESTS_2( bloc, cont, grid )
DEF_MEXPR_TESTS_2( bloc, bloc, cont )
DEF_MEXPR_TESTS_2( bloc, bloc, bloc )
DEF_MEXPR_TESTS_2( bloc, bloc, grid )
DEF_MEXPR_TESTS_2( bloc, grid, cont )
DEF_MEXPR_TESTS_2( bloc, grid, bloc )
DEF_MEXPR_TESTS_2( bloc, grid, grid )

DEF_MEXPR_TESTS_2( grid, cont, cont )
DEF_MEXPR_TESTS_2( grid, cont, bloc )
DEF_MEXPR_TESTS_2( grid, cont, grid )
DEF_MEXPR_TESTS_2( grid, bloc, cont )
DEF_MEXPR_TESTS_2( grid, bloc, bloc )
DEF_MEXPR_TESTS_2( grid, bloc, grid )
DEF_MEXPR_TESTS_2( grid, grid, cont )
DEF_MEXPR_TESTS_2( grid, grid, bloc )
DEF_MEXPR_TESTS_2( grid, grid, grid )

// Ternary expression


#define DEF_MEXPR_TESTS_3( stag1, stag2, stag3, dtag ) \
		MN_CASE( map_expr, ternary_##stag1##_##stag2##_##stag3##_##dtag ) { test_mapexpr_3<stag1, stag2, stag3, dtag, M, N>(); } \
		BEGIN_TPACK( ternary_map_expr_##stag1##_##stag2##_##stag3##_##dtag ) \
		ADD_MN_CASE_3X3( map_expr, ternary_##stag1##_##stag2##_##stag3##_##dtag, DM, DN ) \
		END_TPACK

DEF_MEXPR_TESTS_3( cont, cont, cont, cont )
DEF_MEXPR_TESTS_3( cont, cont, cont, bloc )
DEF_MEXPR_TESTS_3( cont, cont, cont, grid )
DEF_MEXPR_TESTS_3( cont, cont, bloc, cont )
DEF_MEXPR_TESTS_3( cont, cont, bloc, bloc )
DEF_MEXPR_TESTS_3( cont, cont, bloc, grid )
DEF_MEXPR_TESTS_3( cont, cont, grid, cont )
DEF_MEXPR_TESTS_3( cont, cont, grid, bloc )
DEF_MEXPR_TESTS_3( cont, cont, grid, grid )

DEF_MEXPR_TESTS_3( cont, bloc, cont, cont )
DEF_MEXPR_TESTS_3( cont, bloc, cont, bloc )
DEF_MEXPR_TESTS_3( cont, bloc, cont, grid )
DEF_MEXPR_TESTS_3( bloc, cont, bloc, cont )
DEF_MEXPR_TESTS_3( bloc, cont, bloc, bloc )
DEF_MEXPR_TESTS_3( bloc, cont, bloc, grid )
DEF_MEXPR_TESTS_3( cont, bloc, grid, cont )
DEF_MEXPR_TESTS_3( cont, bloc, grid, bloc )
DEF_MEXPR_TESTS_3( cont, bloc, grid, grid )


// Compound expressions

MN_CASE( map_expr, compound_expr )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef mat_host<bloc, double, M, N> shost1_t;
	typedef mat_host<bloc, double, M, N> shost2_t;
	typedef mat_host<cont, double, M, N> shost3_t;
	typedef mat_host<cont, double, M, N> dhost_t;

	typedef typename shost1_t::cmat_t smat1_t;
	typedef typename shost2_t::cmat_t smat2_t;
	typedef typename shost3_t::cmat_t smat3_t;
	typedef typename dhost_t::mat_t dmat_t;

	shost1_t s1_h(m, n);
	shost2_t s2_h(m, n);
	shost3_t s3_h(m, n);
	dhost_t d_h(m, n);

	s1_h.fill_rand();
	s2_h.fill_rand();
	s2_h.fill_rand();

	smat1_t s1 = s1_h.get_cmat();
	smat2_t s2 = s2_h.get_cmat();
	smat3_t s3 = s3_h.get_cmat();
	dmat_t d = d_h.get_mat();

	double cv = 2.5;

	auto e = make_map_expr(
			fma_(),
			s1,
			make_map_expr_fix2(sub_(), s2, cv),
			make_map_expr(sqr_(), s3) );

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n);

	d = e;

	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = s1(i, j) * (s2(i, j) - cv) + s3(i, j) * s3(i, j);
	}

	double tol = 1.0e-14;
	ASSERT_MAT_APPROX(m, n, d, r, tol);
}

BEGIN_TPACK( compound_map_expr )
	ADD_MN_CASE_3X3( map_expr, compound_expr, DM, DN )
END_TPACK



BEGIN_MAIN_SUITE

	// unary

	ADD_TPACK( unary_map_expr_cont_cont )
	ADD_TPACK( unary_map_expr_cont_bloc )
	ADD_TPACK( unary_map_expr_cont_grid )
	ADD_TPACK( unary_map_expr_bloc_cont )
	ADD_TPACK( unary_map_expr_bloc_bloc )
	ADD_TPACK( unary_map_expr_bloc_grid )
	ADD_TPACK( unary_map_expr_grid_cont )
	ADD_TPACK( unary_map_expr_grid_bloc )
	ADD_TPACK( unary_map_expr_grid_grid )

	// binary

	ADD_TPACK( binary_map_expr_cont_cont_cont )
	ADD_TPACK( binary_map_expr_cont_cont_bloc )
	ADD_TPACK( binary_map_expr_cont_cont_grid )
	ADD_TPACK( binary_map_expr_cont_bloc_cont )
	ADD_TPACK( binary_map_expr_cont_bloc_bloc )
	ADD_TPACK( binary_map_expr_cont_bloc_grid )
	ADD_TPACK( binary_map_expr_cont_grid_cont )
	ADD_TPACK( binary_map_expr_cont_grid_bloc )
	ADD_TPACK( binary_map_expr_cont_grid_grid )

	ADD_TPACK( binary_map_expr_bloc_cont_cont )
	ADD_TPACK( binary_map_expr_bloc_cont_bloc )
	ADD_TPACK( binary_map_expr_bloc_cont_grid )
	ADD_TPACK( binary_map_expr_bloc_bloc_cont )
	ADD_TPACK( binary_map_expr_bloc_bloc_bloc )
	ADD_TPACK( binary_map_expr_bloc_bloc_grid )
	ADD_TPACK( binary_map_expr_bloc_grid_cont )
	ADD_TPACK( binary_map_expr_bloc_grid_bloc )
	ADD_TPACK( binary_map_expr_bloc_grid_grid )

	ADD_TPACK( binary_map_expr_grid_cont_cont )
	ADD_TPACK( binary_map_expr_grid_cont_bloc )
	ADD_TPACK( binary_map_expr_grid_cont_grid )
	ADD_TPACK( binary_map_expr_grid_bloc_cont )
	ADD_TPACK( binary_map_expr_grid_bloc_bloc )
	ADD_TPACK( binary_map_expr_grid_bloc_grid )
	ADD_TPACK( binary_map_expr_grid_grid_cont )
	ADD_TPACK( binary_map_expr_grid_grid_bloc )
	ADD_TPACK( binary_map_expr_grid_grid_grid )

	// ternary

	ADD_TPACK( ternary_map_expr_cont_cont_cont_cont )
	ADD_TPACK( ternary_map_expr_cont_cont_cont_bloc )
	ADD_TPACK( ternary_map_expr_cont_cont_cont_grid )
	ADD_TPACK( ternary_map_expr_cont_cont_bloc_cont )
	ADD_TPACK( ternary_map_expr_cont_cont_bloc_bloc )
	ADD_TPACK( ternary_map_expr_cont_cont_bloc_grid )
	ADD_TPACK( ternary_map_expr_cont_cont_grid_cont )
	ADD_TPACK( ternary_map_expr_cont_cont_grid_bloc )
	ADD_TPACK( ternary_map_expr_cont_cont_grid_grid )

	ADD_TPACK( ternary_map_expr_cont_bloc_cont_cont )
	ADD_TPACK( ternary_map_expr_cont_bloc_cont_bloc )
	ADD_TPACK( ternary_map_expr_cont_bloc_cont_grid )
	ADD_TPACK( ternary_map_expr_bloc_cont_bloc_cont )
	ADD_TPACK( ternary_map_expr_bloc_cont_bloc_bloc )
	ADD_TPACK( ternary_map_expr_bloc_cont_bloc_grid )
	ADD_TPACK( ternary_map_expr_cont_bloc_grid_cont )
	ADD_TPACK( ternary_map_expr_cont_bloc_grid_bloc )
	ADD_TPACK( ternary_map_expr_cont_bloc_grid_grid )

	// compound

	ADD_TPACK( compound_map_expr )

END_MAIN_SUITE

