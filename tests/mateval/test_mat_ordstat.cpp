/**
 * @file test_mat_ordstat.cpp
 *
 * @brief Unit testing for Order statistics
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/matrix_sort.h>
#include <light_mat/mateval/matrix_ordstats.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 15;
const index_t DM2 = 12;
const index_t DN = 16;

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& a)
{
	for (index_t i = 0; i < a.nelems(); ++i)
	{
		a[i] = (double)std::rand() / (double)RAND_MAX;
	}
}


SIMPLE_CASE( vec_find_max_min )
{
	const index_t n = DM;

	dense_col<double> a(n);
	fill_ran(a);

	dense_col<index_t> si = sorted_idx(a);

	index_t imin0 = si[0];
	index_t imax0 = si[n-1];
	double vmin0 = a[imin0];
	double vmax0 = a[imax0];

	index_t imax = find_imax(a);
	index_t imin = find_imin(a);

	ASSERT_EQ( imax, imax0 );
	ASSERT_EQ( imin, imin0 );

	imax = -1;
	imin = -1;
	double vmax = 0;
	double vmin = 0;

	std::tie( imax, vmax ) = find_max(a);
	std::tie( imin, vmin ) = find_min(a);

	ASSERT_EQ( imax, imax0 );
	ASSERT_EQ( imin, imin0 );
	ASSERT_EQ( vmax, vmax0 );
	ASSERT_EQ( vmin, vmin0 );
}


SIMPLE_CASE( colwise_find_max_min )
{
	const index_t m = DM;
	const index_t n = DN;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	dense_matrix<index_t> si = colwise_sorted_idx(a);
	dense_matrix<double> sx = colwise_sorted(a);

	dense_row<index_t>  ri_max(n, zero());
	dense_row<double> rx_max(n, zero());
	dense_row<index_t>  ri_min(n, zero());
	dense_row<double> rx_min(n, zero());

	colwise_find_imax(a, ri_max);
	colwise_find_imin(a, ri_min);

	ASSERT_VEC_EQ( n, ri_max, si.row(m - 1) );
	ASSERT_VEC_EQ( n, ri_min, si.row(0) );

	zero(ri_max);
	zero(ri_min);

	colwise_find_max(a, ri_max, rx_max );
	colwise_find_min(a, ri_min, rx_min );

	ASSERT_VEC_EQ( n, ri_max, si.row(m - 1) );
	ASSERT_VEC_EQ( n, ri_min, si.row(0) );
	ASSERT_VEC_EQ( n, rx_max, sx.row(m - 1) );
	ASSERT_VEC_EQ( n, rx_min, sx.row(0) );
}


SIMPLE_CASE( vec_nth_elem )
{
	const index_t n = DM;

	dense_col<double> a(n);
	fill_ran(a);

	dense_col<double> sx = sorted(a);

	for (index_t k = 0; k < n; ++k)
	{
		double r = nth_element(a, k);
		ASSERT_EQ( r, sx[k] );
	}
}

SIMPLE_CASE( colwise_nth_elem )
{
	const index_t m = DM;
	const index_t n = DN;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	dense_matrix<double> sx = colwise_sorted(a);

	for (index_t k = 0; k < m; ++k)
	{
		dense_row<double> r(n, zero());
		colwise_nth_element(a, k, r);
		ASSERT_VEC_EQ( n, r, sx.row(k) );
	}
}


SIMPLE_CASE( vec_median_odd )
{
	const index_t n = DM;
	const index_t mid = (n - 1) / 2;

	dense_col<double> a(n);
	fill_ran(a);

	dense_col<double> sx = sorted(a);

	double r = median(a);
	ASSERT_EQ( r, sx[mid] );
}

SIMPLE_CASE( vec_median_even )
{
	const index_t n = DM2;
	const index_t mid1 = n / 2;
	const index_t mid0 = mid1 - 1;

	dense_col<double> a(n);
	fill_ran(a);

	dense_col<double> sx = sorted(a);
	double r0 = (sx[mid0] + sx[mid1]) / 2;

	double r = median(a);
	ASSERT_APPROX( r, r0, 1.0e-15 );
}

SIMPLE_CASE( colwise_median_odd )
{
	const index_t m = DM;
	const index_t n = DN;
	const index_t mid = (m - 1) / 2;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	dense_matrix<double> sx = colwise_sorted(a);

	dense_row<double> r(n, zero());
	colwise_median(a, r);
	ASSERT_VEC_EQ( n, r, sx.row(mid) );
}

SIMPLE_CASE( colwise_median_even )
{
	const index_t m = DM2;
	const index_t n = DN;
	const index_t mid1 = m / 2;
	const index_t mid0 = mid1 - 1;

	dense_matrix<double> a(m, n);
	fill_ran(a);

	dense_matrix<double> sx = colwise_sorted(a);
	dense_row<double> r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		r0[j] = (sx(mid0, j) + sx(mid1, j)) / 2;
	}

	dense_row<double> r(n, zero());
	colwise_median(a, r);

	ASSERT_VEC_APPROX( n, r, r0, 1.0e-15 );
}


AUTO_TPACK( test_find_max_min )
{
	ADD_SIMPLE_CASE( vec_find_max_min )
	ADD_SIMPLE_CASE( colwise_find_max_min )
}

AUTO_TPACK( test_nth_elem )
{
	ADD_SIMPLE_CASE( vec_nth_elem )
	ADD_SIMPLE_CASE( colwise_nth_elem )
}

AUTO_TPACK( test_median )
{
	ADD_SIMPLE_CASE( vec_median_odd )
	ADD_SIMPLE_CASE( vec_median_even )
	ADD_SIMPLE_CASE( colwise_median_odd )
	ADD_SIMPLE_CASE( colwise_median_even )
}




