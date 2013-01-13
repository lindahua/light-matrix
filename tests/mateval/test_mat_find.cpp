/**
 * @file test_mat_find.cpp
 *
 * @brief Unit testing of matrix finding
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_pred.h>
#include <light_mat/mateval/matrix_find.h>
#include <vector>
#include <algorithm>

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


MN_CASE( mat_count )
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

MN_CASE( mat_count_ex )
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


MN_CASE( mat_colwise_count )
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


MN_CASE( mat_colwise_count_ex )
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


SIMPLE_CASE( mat_findl )
{
	const index_t m = DM;
	const index_t n = DN;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	double t = 0.4;

	std::vector<index_t> vecl;
	findl_to(a > t, vecl);

	std::vector<index_t> vecl0;
	for (index_t i = 0; i < m * n; ++i)
	{
		if (a[i] > t) vecl0.push_back(i);
	}

	ASSERT_EQ( vecl.size(), vecl0.size() );
	ASSERT_TRUE( std::equal(vecl.begin(), vecl.end(), vecl0.begin()) );
}


SIMPLE_CASE( mat_findij )
{
	const index_t m = DM;
	const index_t n = DN;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	double t = 0.4;

	std::vector<index_t> veci;
	std::vector<index_t> vecj;
	find_to(a > t, veci, vecj);

	std::vector<index_t> veci0;
	std::vector<index_t> vecj0;

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (a(i, j) > t)
			{
				veci0.push_back(i);
				vecj0.push_back(j);
			}
		}
	}

	ASSERT_EQ( veci.size(), veci0.size() );
	ASSERT_EQ( vecj.size(), vecj0.size() );
	ASSERT_TRUE( std::equal(veci.begin(), veci.end(), veci0.begin()) );
	ASSERT_TRUE( std::equal(vecj.begin(), vecj.end(), vecj0.begin()) );
}



AUTO_TPACK( mat_count )
{
	ADD_MN_CASE_3X3( mat_count, DM, DN )
}

AUTO_TPACK( mat_count_ex )
{
	ADD_MN_CASE_3X3( mat_count_ex, DM, DN )
}

AUTO_TPACK( mat_colwise_count )
{
	ADD_MN_CASE_3X3( mat_colwise_count, DM, DN )
}

AUTO_TPACK( mat_colwise_count_ex )
{
	ADD_MN_CASE_3X3( mat_colwise_count_ex, DM, DN )
}

AUTO_TPACK( mat_find )
{
	ADD_SIMPLE_CASE( mat_findl )
	ADD_SIMPLE_CASE( mat_findij )
}





