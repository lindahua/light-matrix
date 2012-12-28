/**
 * @file test_more_reduce.cpp
 *
 * @brief Unit testing of more reduction functions
 *
 * @author Dahua Lin
 */


#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_reduce.h>

#include <light_mat/mateval/mat_minmax.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;


inline double randunif()
{
	double u = (double)std::rand() / double(RAND_MAX);
	return u * 2.0 - 1.0;
}

template<class Mat, typename T>
void fill_rand(IRegularMatrix<Mat, T>& mat)
{
	for (index_t j = 0; j < mat.ncolumns(); ++j)
	{
		for (index_t i = 0; i < mat.nrows(); ++i)
		{
			mat(i, j) = randunif();
		}
	}
}

const index_t max_len = 48;
const index_t max_nrows = 32;
index_t test_nrows[] = { 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 28, 32 };
const unsigned int ntest_nrows = sizeof(test_nrows) / sizeof(index_t);

SIMPLE_CASE( minmax, full )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		double rmin = minimum(s(range(0, k)));
		double rmax = maximum(s(range(0, k)));

		minmax_stat<double> rs = minmax(s(range(0, k)));

		ASSERT_EQ(rs.min_value, rmin);
		ASSERT_EQ(rs.max_value, rmax);
	}
}

SIMPLE_CASE( minmax, colwise )
{
	const index_t n = 6;
	dense_matrix<double> src(max_nrows, n);
	fill_rand(src);

	dense_row<double> min_r0(n, zero());
	dense_row<double> max_r0(n, zero());
	dense_row<double> min_r(n, zero());
	dense_row<double> max_r(n, zero());

	for (unsigned k = 0; k < ntest_nrows; ++k)
	{
		index_t cl = test_nrows[k];
		ref_block<double> s = src(range(0, cl), whole());

		colwise_minimum(s, min_r0);
		colwise_maximum(s, max_r0);
		colwise_minmax(s, min_r, max_r);

		ASSERT_MAT_EQ(1, n, min_r, min_r0);
		ASSERT_MAT_EQ(1, n, max_r, max_r0);
	}
}



BEGIN_TPACK( mat_minmax )
	ADD_SIMPLE_CASE( minmax, full )
	ADD_SIMPLE_CASE( minmax, colwise )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_minmax )
END_MAIN_SUITE


