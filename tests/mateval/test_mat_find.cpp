/**
 * @file test_mat_find.cpp
 *
 * @brief Unit testing of matrix finding
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_pred.h>
#include <light_mat/mateval/matrix_find.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 9;
const index_t DN = 8;
const index_t LDim = 12;

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& a)
{
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		a[i] = (double)std::rand() / (double)RAND_MAX;
	}
}


MN_CASE( mat_count, count )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	double t = 0.4;
	size_t r = count(a > t);

	size_t r0 = 0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			if (a(i, j) > t) ++r0;
	}

	ASSERT_EQ( r, r0 );
}

MN_CASE( mat_count, count_ex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double> a0(LDim, n);
	fill_ran(a0);

	ref_block<double, M, N> a(a0.ptr_data(), m, n, LDim);

	double t = 0.4;
	size_t r = count(a > t);

	size_t r0 = 0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			if (a(i, j) > t) ++r0;
	}

	ASSERT_EQ( r, r0 );
}


MN_CASE( mat_count, colwise_count )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	fill_ran(a);

	double t = 0.4;
	dense_row<index_t, N> r(n);
	colwise_count(a > t, r);

	dense_row<index_t, N> r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		index_t rj = 0;
		for (index_t i = 0; i < m; ++i)
			if (a(i, j) > t) ++rj;
		r0[j] = rj;
	}

	ASSERT_VEC_EQ( n, r, r0 );
}


MN_CASE( mat_count, colwise_count_ex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<double> a0(LDim, n);
	fill_ran(a0);
	ref_block<double, M, N> a(a0.ptr_data(), m, n, LDim);

	double t = 0.4;
	dense_row<index_t, N> r(n);
	colwise_count(a > t, r);

	dense_row<index_t, N> r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		index_t rj = 0;
		for (index_t i = 0; i < m; ++i)
			if (a(i, j) > t) ++rj;
		r0[j] = rj;
	}

	ASSERT_VEC_EQ( n, r, r0 );
}




BEGIN_TPACK( mat_count )
	ADD_MN_CASE_3X3( mat_count, count, DM, DN )
END_TPACK

BEGIN_TPACK( mat_count_ex )
	ADD_MN_CASE_3X3( mat_count, count_ex, DM, DN )
END_TPACK

BEGIN_TPACK( mat_colwise_count )
	ADD_MN_CASE_3X3( mat_count, colwise_count, DM, DN )
END_TPACK

BEGIN_TPACK( mat_colwise_count_ex )
	ADD_MN_CASE_3X3( mat_count, colwise_count_ex, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_count )
	ADD_TPACK( mat_count_ex )
	ADD_TPACK( mat_colwise_count )
	ADD_TPACK( mat_colwise_count_ex )
END_MAIN_SUITE






