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

	typedef map_expr<sqr_, smat1_t> expr_t;
	expr_t e = make_map_expr(sqr_(), s1);

	ASSERT_EQ( e.nrows(), m );
	ASSERT_EQ( e.ncolumns(), n );
	ASSERT_EQ( e.nelems(), m * n);

	d = e;

	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			r(i, j) = s1(i, j) + s2(i, j);
	}

	double tol = 1.0e-14;
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

MN_CASE( map_expr, binary_cont_cont_cont ) { test_mapexpr_2<cont, cont, cont, M, N>(); }


BEGIN_MAIN_SUITE
	ADD_TPACK( unary_map_expr_cont_cont )
	ADD_TPACK( unary_map_expr_cont_bloc )
	ADD_TPACK( unary_map_expr_cont_grid )
	ADD_TPACK( unary_map_expr_bloc_cont )
	ADD_TPACK( unary_map_expr_bloc_bloc )
	ADD_TPACK( unary_map_expr_bloc_grid )
	ADD_TPACK( unary_map_expr_grid_cont )
	ADD_TPACK( unary_map_expr_grid_bloc )
	ADD_TPACK( unary_map_expr_grid_grid )
END_MAIN_SUITE

