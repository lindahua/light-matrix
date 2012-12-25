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

	ASSERT_MAT_EQ(m, n, d, r);
}


// Unary expressions

MN_CASE( map_expr, unary_cont_cont ) { test_mapexpr_1<cont, cont, M, N>(); }
MN_CASE( map_expr, unary_cont_bloc ) { test_mapexpr_1<cont, bloc, M, N>(); }
MN_CASE( map_expr, unary_cont_grid ) { test_mapexpr_1<cont, grid, M, N>(); }
MN_CASE( map_expr, unary_bloc_cont ) { test_mapexpr_1<bloc, cont, M, N>(); }
MN_CASE( map_expr, unary_bloc_bloc ) { test_mapexpr_1<bloc, bloc, M, N>(); }
MN_CASE( map_expr, unary_bloc_grid ) { test_mapexpr_1<bloc, grid, M, N>(); }
MN_CASE( map_expr, unary_grid_cont ) { test_mapexpr_1<grid, cont, M, N>(); }
MN_CASE( map_expr, unary_grid_bloc ) { test_mapexpr_1<grid, bloc, M, N>(); }
MN_CASE( map_expr, unary_grid_grid ) { test_mapexpr_1<grid, grid, M, N>(); }

BEGIN_TPACK( unary_map_expr_cont_cont )
ADD_MN_CASE_3X3( map_expr, unary_cont_cont, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_cont_bloc )
ADD_MN_CASE_3X3( map_expr, unary_cont_bloc, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_cont_grid )
ADD_MN_CASE_3X3( map_expr, unary_cont_grid, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_bloc_cont )
ADD_MN_CASE_3X3( map_expr, unary_bloc_cont, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_bloc_bloc )
ADD_MN_CASE_3X3( map_expr, unary_bloc_bloc, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_bloc_grid )
ADD_MN_CASE_3X3( map_expr, unary_bloc_grid, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_grid_cont )
ADD_MN_CASE_3X3( map_expr, unary_grid_cont, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_grid_bloc )
ADD_MN_CASE_3X3( map_expr, unary_grid_bloc, DM, DN )
END_TPACK
BEGIN_TPACK( unary_map_expr_grid_grid )
ADD_MN_CASE_3X3( map_expr, unary_grid_grid, DM, DN )
END_TPACK


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

