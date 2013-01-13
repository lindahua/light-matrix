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


template<class STag, int M, int N>
void test_asrow_view()
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

	dense_matrix<double> r2(1, m * n);
	do_fill_lin(r2.ptr_data(), m * n);

	const cmat_t& Ac = A;
	const mat_t& Bc = B;

	auto ac = as_row(Ac);
	auto bc = as_row(Bc);

	auto a = as_row(A);
	auto b = as_row(B);

	ASSERT_EQ( ac.nrows(), 1 );
	ASSERT_EQ( ac.ncolumns(), m * n );

	ASSERT_EQ( bc.nrows(), 1 );
	ASSERT_EQ( bc.ncolumns(), m * n );

	ASSERT_EQ( a.nrows(), 1 );
	ASSERT_EQ( a.ncolumns(), m * n );

	ASSERT_EQ( b.nrows(), 1 );
	ASSERT_EQ( b.ncolumns(), m * n );

	ASSERT_VEC_EQ( m * n, ac, r );
	ASSERT_VEC_EQ( m * n, bc, r );
	ASSERT_VEC_EQ( m * n, a, r );
	ASSERT_VEC_EQ( m * n, b, r );

	as_row(B) = r2;
	ASSERT_VEC_EQ( m * n, b, r2 );
}


SIMPLE_CASE( vector_as_col )
{
	std::vector<double> x;
	const index_t n = 10;
	dense_col<double> r(n);

	for (index_t i = 0; i < 10; ++i)
	{
		r[i] = double(i+1);
		x.push_back(r[i]);
	}

	auto xv = as_col(x);

	ASSERT_EQ( xv.nrows(), n );
	ASSERT_EQ( xv.ncolumns(), 1 );
	ASSERT_EQ( xv.nelems(), n );

	ASSERT_VEC_EQ( n, xv, r );
}

SIMPLE_CASE( vector_as_row )
{
	std::vector<double> x;
	const index_t n = 10;
	dense_row<double> r(n);

	for (index_t i = 0; i < 10; ++i)
	{
		r[i] = double(i+1);
		x.push_back(r[i]);
	}

	auto xv = as_row(x);

	ASSERT_EQ( xv.nrows(), 1 );
	ASSERT_EQ( xv.ncolumns(), n );
	ASSERT_EQ( xv.nelems(), n );

	ASSERT_VEC_EQ( n, xv, r );
}



MN_CASE( mat_as_col_cont )
{
	test_ascol_view<cont, M, N>();
}

MN_CASE( mat_as_col_bloc )
{
	test_ascol_view<bloc, M, N>();
}

MN_CASE( mat_as_col_grid )
{
	test_ascol_view<grid, M, N>();
}

MN_CASE( mat_as_row_cont )
{
	test_asrow_view<cont, M, N>();
}

MN_CASE( mat_as_row_bloc )
{
	test_asrow_view<bloc, M, N>();
}

MN_CASE( mat_as_row_grid )
{
	test_asrow_view<grid, M, N>();
}


LTEST_INIT_AUTOSUITE

AUTO_TPACK( vector_adapter )
{
	ADD_SIMPLE_CASE( vector_as_col )
	ADD_SIMPLE_CASE( vector_as_row )
}

AUTO_TPACK( mat_as_col )
{
	ADD_MN_CASE_3X3( mat_as_col_cont, DM, DN )

	ADD_MN_CASE( mat_as_col_bloc, 1, 1 )
	ADD_MN_CASE( mat_as_col_bloc, 0, 1 )
	ADD_MN_CASE( mat_as_col_bloc, DM, 1 )
	ADD_MN_CASE( mat_as_col_bloc, 1, 0 )
	ADD_MN_CASE( mat_as_col_bloc, 1, DN )

	ADD_MN_CASE( mat_as_col_grid, 1, 1 )
	ADD_MN_CASE( mat_as_col_grid, 0, 1 )
	ADD_MN_CASE( mat_as_col_grid, DM, 1 )
	ADD_MN_CASE( mat_as_col_grid, 1, 0 )
	ADD_MN_CASE( mat_as_col_grid, 1, DN )
}


AUTO_TPACK( mat_as_row )
{
	ADD_MN_CASE_3X3( mat_as_row_cont, DM, DN )

	ADD_MN_CASE( mat_as_row_bloc, 1, 1 )
	ADD_MN_CASE( mat_as_row_bloc, 0, 1 )
	ADD_MN_CASE( mat_as_row_bloc, DM, 1 )
	ADD_MN_CASE( mat_as_row_bloc, 1, 0 )
	ADD_MN_CASE( mat_as_row_bloc, 1, DN )

	ADD_MN_CASE( mat_as_row_grid, 1, 1 )
	ADD_MN_CASE( mat_as_row_grid, 0, 1 )
	ADD_MN_CASE( mat_as_row_grid, DM, 1 )
	ADD_MN_CASE( mat_as_row_grid, 1, 0 )
	ADD_MN_CASE( mat_as_row_grid, 1, DN )
}





