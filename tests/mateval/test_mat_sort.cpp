/**
 * @file test_mat_sort.cpp
 *
 * @brief Testing of matrix sorting
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/matrix_sort.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 15;
const index_t DN = 8;

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& a)
{
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		a[i] = (double)std::rand() / (double)RAND_MAX;
	}
}


// core testing function

template<class A, class R>
bool test_sorted(const A& a, const R& r, asc_)
{
	if (!have_same_shape(a, r)) return false;

	index_t m = a.nrows();
	index_t n = a.ncolumns();
	dense_matrix<double> c(a);

	std::sort(begin(c), end(c), std::less<double>());
	return is_sorted(r, asc_()) && ltest::test_matrix_equal(m, n, r, c);
}

template<class A, class R>
bool test_sorted(const A& a, const R& r, desc_)
{
	if (!have_same_shape(a, r)) return false;

	index_t m = a.nrows();
	index_t n = a.ncolumns();
	dense_matrix<double> c(a);

	std::sort(begin(c), end(c), std::greater<double>());
	return is_sorted(r, desc_()) && ltest::test_matrix_equal(m, n, r, c);
}

template<class A, class R>
bool test_cw_sorted(const A& a, const R& r, asc_)
{
	if (!have_same_shape(a, r)) return false;

	index_t n = a.ncolumns();
	for (index_t j = 0; j < n; ++j)
	{
		if (!test_sorted(a.column(j), r.column(j), asc_()))
			return false;
	}

	return true;
}

template<class A, class R>
bool test_cw_sorted(const A& a, const R& r, desc_)
{
	if (!have_same_shape(a, r)) return false;

	index_t n = a.ncolumns();
	for (index_t j = 0; j < n; ++j)
	{
		if (!test_sorted(a.column(j), r.column(j), desc_()))
			return false;
	}

	return true;
}


template<class A, class R, typename S>
bool test_sorted_idx(const A& a, const R& ri, S)
{
	if (!have_same_shape(a, ri)) return false;

	dense_matrix<double> r(a.nrows(), a.ncolumns());
	for (index_t i = 0; i < a.nelems(); ++i) r[i] = a[ri[i]];

	return test_sorted(a, r, S());
}

template<class A, class R, typename S>
bool test_cw_sorted_idx(const A& a, const R& ri, S)
{
	if (!have_same_shape(a, ri)) return false;

	index_t m = a.nrows();
	index_t n = a.ncolumns();
	dense_matrix<double> r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) r(i, j) = a(ri(i, j), j);

	return test_cw_sorted(a, r, S());
}


template<class A, class R, typename S>
bool test_sorted_ex(const A& a, const R& re, S)
{
	if (!have_same_shape(a, re)) return false;

	dense_matrix<double> rv;
	dense_matrix<index_t> ri;
	unzip(re, rv, ri);

	return test_sorted(a, rv, S()) && test_sorted_idx(a, ri, S());
}

template<class A, class R, typename S>
bool test_cw_sorted_ex(const A& a, const R& re, S)
{
	if (!have_same_shape(a, re)) return false;

	dense_matrix<double> rv;
	dense_matrix<index_t> ri;
	unzip(re, rv, ri);

	return test_cw_sorted(a, rv, S()) && test_cw_sorted_idx(a, ri, S());
}



// Test cases

MN_CASE( mat_inplace_sort )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	sort(a, asc_());
	ASSERT_TRUE( test_sorted(a0, a, asc_()) );

	fill_ran(a);
	copy(a, a0);
	sort(a, desc_());
	ASSERT_TRUE( test_sorted(a0, a, desc_()) );

	fill_ran(a);
	copy(a, a0);
	sort(a);
	ASSERT_TRUE( test_sorted(a0, a, asc_()) );
}

MN_CASE( mat_colwise_inplace_sort )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	colwise_sort(a, asc_());
	ASSERT_TRUE( test_cw_sorted(a0, a, asc_()) );

	fill_ran(a);
	copy(a, a0);
	colwise_sort(a, desc_());
	ASSERT_TRUE( test_cw_sorted(a0, a, desc_()) );

	fill_ran(a);
	copy(a, a0);
	colwise_sort(a);
	ASSERT_TRUE( test_cw_sorted(a0, a, asc_()) );
}

MN_CASE( mat_copy_sort )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b1 = sorted(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b2 = sorted(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b3 = sorted(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted(a, b3, asc_()) );
}

MN_CASE( mat_colwise_copy_sort )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b1 = colwise_sorted(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );

	ASSERT_TRUE( test_cw_sorted(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b2 = colwise_sorted(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<double, M, N> b3 = colwise_sorted(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted(a, b3, asc_()) );
}


MN_CASE( mat_sort_idx )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b1 = sorted_idx(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_idx(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b2 = sorted_idx(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_idx(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b3 = sorted_idx(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_idx(a, b3, asc_()) );
}

MN_CASE( mat_colwise_sort_idx )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b1 = colwise_sorted_idx(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_idx(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b2 = colwise_sorted_idx(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_idx(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<index_t, M, N> b3 = colwise_sorted_idx(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_idx(a, b3, asc_()) );
}


MN_CASE( mat_sort_ex )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());
	typedef std::pair<double, index_t> xt;

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b1 = sorted_ex(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_ex(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b2 = sorted_ex(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_ex(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b3 = sorted_ex(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_sorted_ex(a, b3, asc_()) );
}


MN_CASE( mat_colwise_sort_ex )
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	dense_matrix<double, M, N> a(m, n);
	dense_matrix<double, M, N> a0(m, n, zero());
	typedef std::pair<double, index_t> xt;

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b1 = colwise_sorted_ex(a, asc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_ex(a, b1, asc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b2 = colwise_sorted_ex(a, desc_());
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_ex(a, b2, desc_()) );

	fill_ran(a);
	copy(a, a0);
	dense_matrix<xt, M, N> b3 = colwise_sorted_ex(a);
	ASSERT_MAT_EQ( m, n, a, a0 );
	ASSERT_TRUE( test_cw_sorted_ex(a, b3, asc_()) );
}



LTEST_INIT_AUTOSUITE

AUTO_TPACK( mat_inplace_sort )
{
	ADD_MN_CASE_3X3( mat_inplace_sort, DM, DN )
}

AUTO_TPACK( mat_colwise_inplace_sort )
{
	ADD_MN_CASE_3X3( mat_colwise_inplace_sort, DM, DN )
}

AUTO_TPACK( mat_copy_sort )
{
	ADD_MN_CASE_3X3( mat_copy_sort, DM, DN )
}

AUTO_TPACK( mat_colwise_copy_sort )
{
	ADD_MN_CASE_3X3( mat_colwise_copy_sort, DM, DN )
}

AUTO_TPACK( mat_sort_idx )
{
	ADD_MN_CASE_3X3( mat_sort_idx, DM, DN )
}

AUTO_TPACK( mat_colwise_sort_idx )
{
	ADD_MN_CASE_3X3( mat_colwise_sort_idx, DM, DN )
}

AUTO_TPACK( mat_sort_ex )
{
	ADD_MN_CASE_3X3( mat_sort_ex, DM, DN )
}

AUTO_TPACK( mat_colwise_sort_ex )
{
	ADD_MN_CASE_3X3( mat_colwise_sort_ex, DM, DN )
}





