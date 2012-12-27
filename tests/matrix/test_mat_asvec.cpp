/**
 * @file test_mat_asvec.cpp
 *
 * @brief Unit testing of matrix as-vector views
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_asvec.h>

using namespace lmat;
using namespace lmat::test;


template<class STag, int M, int N>
void test_ascol_view()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_host<STag, double, M, N> host_t;
	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	host_t host(m, n);
	host.fill_rand();

	cmat_t A = host.get_cmat();
	mat_t B = host.get_mat();
	dense_matrix<double> r = A;

	dense_matrix<double> r2(m * n, 1);
	do_fill_lin(r2.ptr_data(), m * n);

	const cmat_t& Ac = A;
	const mat_t& Bc = B;

	auto ac = as_col(Ac);
	auto bc = as_col(Bc);

	auto a = as_col(A);
	auto b = as_col(B);

	ASSERT_EQ( ac.nrows(), m * n );
	ASSERT_EQ( ac.ncolumns(), 1 );

	ASSERT_EQ( bc.nrows(), m * n );
	ASSERT_EQ( bc.ncolumns(), 1 );

	ASSERT_EQ( a.nrows(), m * n );
	ASSERT_EQ( a.ncolumns(), 1 );

	ASSERT_EQ( b.nrows(), m * n );
	ASSERT_EQ( b.ncolumns(), 1 );

	ASSERT_VEC_EQ( m * n, ac, r );
	ASSERT_VEC_EQ( m * n, bc, r );
	ASSERT_VEC_EQ( m * n, a, r );
	ASSERT_VEC_EQ( m * n, b, r );

	as_col(B) = r2;
	ASSERT_VEC_EQ( m * n, b, r2 );
}


#define DEF_ASVEC_TESTS( V, STag ) \
	MN_CASE( mat_as_##V, STag##_as_##V ) { test_as##V##_view<STag, M, N>(); }

DEF_ASVEC_TESTS( col, cont )
DEF_ASVEC_TESTS( col, bloc )
DEF_ASVEC_TESTS( col, grid )


BEGIN_TPACK( mat_as_col_cont )
	ADD_MN_CASE_3X3( mat_as_col, cont_as_col, DM, DN )
END_TPACK

BEGIN_TPACK( mat_as_col_bloc )
	ADD_MN_CASE( mat_as_col, bloc_as_col, 1, 1 )
	ADD_MN_CASE( mat_as_col, bloc_as_col, 0, 1 )
	ADD_MN_CASE( mat_as_col, bloc_as_col, DM, 1 )
	ADD_MN_CASE( mat_as_col, bloc_as_col, 1, 0 )
	ADD_MN_CASE( mat_as_col, bloc_as_col, 1, DN )
END_TPACK

BEGIN_TPACK( mat_as_col_grid )
	ADD_MN_CASE( mat_as_col, grid_as_col, 1, 1 )
	ADD_MN_CASE( mat_as_col, grid_as_col, 0, 1 )
	ADD_MN_CASE( mat_as_col, grid_as_col, DM, 1 )
	ADD_MN_CASE( mat_as_col, grid_as_col, 1, 0 )
	ADD_MN_CASE( mat_as_col, grid_as_col, 1, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_as_col_cont )
	ADD_TPACK( mat_as_col_bloc )
	ADD_TPACK( mat_as_col_grid )
END_MAIN_SUITE




