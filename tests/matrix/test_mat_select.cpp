/**
 * @file test_mat_select.cpp
 *
 * @brief Unit testing of matrix selection
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

const index_t M0 = 45;
const index_t N0 = 60;

const index_t DM = 9;
const index_t DN = 8;
const index_t LDim = 12;


template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = (double(std::rand()) / RAND_MAX);
	}
}

template<int M, int N>
void fill_randi(dense_matrix<index_t, M, N>& X, index_t U)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = (index_t)std::rand() % U;
	}
}


MN_CASE( mat_select, selectl )
{
	dense_matrix<double> s(M0, N0);
	fill_ran(s);

	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<index_t, M, N> L(m, n);
	fill_randi(L, M0 * N0);

	dense_matrix<double, M, N> r = selectl(s, L);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );

	dense_matrix<double> r0(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r0(i, j) = s[L(i, j)];
	}

	ASSERT_MAT_EQ( m, n, r, r0 );
}

MN_CASE( mat_select, selectl_ex )
{
	dense_matrix<double> s(M0, N0);
	fill_ran(s);

	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<index_t> L0(LDim, n);
	fill_randi(L0, M0 * N0);
	cref_block<index_t, M, N> L(L0.ptr_data(), m, n, LDim);

	dense_matrix<double, M, N> r = selectl(s, L);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );

	dense_matrix<double> r0(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r0(i, j) = s[L(i, j)];
	}

	ASSERT_MAT_EQ( m, n, r, r0 );
}

MN_CASE( mat_select, selectl2 )
{
	dense_matrix<double> s(M0, N0);
	fill_ran(s);

	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<index_t, M, N> I(m, n);
	dense_matrix<index_t, M, N> J(m, n);
	fill_randi(I, M0);
	fill_randi(J, N0);

	dense_matrix<double, M, N> r = selectl(s, I, J);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );

	dense_matrix<double> r0(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r0(i, j) = s(I(i, j), J(i, j));
	}

	ASSERT_MAT_EQ( m, n, r, r0 );
}

MN_CASE( mat_select, selectl2_ex )
{
	dense_matrix<double> s(M0, N0);
	fill_ran(s);

	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<index_t> I0(LDim, n);
	dense_matrix<index_t> J0(LDim, n);
	fill_randi(I0, M0);
	fill_randi(J0, N0);

	cref_block<index_t, M, N> I(I0.ptr_data(), m, n, LDim);
	cref_block<index_t, M, N> J(J0.ptr_data(), m, n, LDim);

	dense_matrix<double, M, N> r = selectl(s, I, J);

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), n );

	dense_matrix<double> r0(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r0(i, j) = s(I(i, j), J(i, j));
	}

	ASSERT_MAT_EQ( m, n, r, r0 );
}



BEGIN_TPACK( mat_selectl )
	ADD_MN_CASE_3X3( mat_select, selectl, DM, DN )
END_TPACK

BEGIN_TPACK( mat_selectl_ex )
	ADD_MN_CASE_3X3( mat_select, selectl_ex, DM, DN )
END_TPACK

BEGIN_TPACK( mat_selectl2 )
	ADD_MN_CASE_3X3( mat_select, selectl2, DM, DN )
END_TPACK

BEGIN_TPACK( mat_selectl2_ex )
	ADD_MN_CASE_3X3( mat_select, selectl2, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_selectl )
	ADD_TPACK( mat_selectl_ex )
	ADD_TPACK( mat_selectl2 )
	ADD_TPACK( mat_selectl2_ex )
END_MAIN_SUITE




