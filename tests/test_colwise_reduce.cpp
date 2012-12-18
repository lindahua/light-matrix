/**
 * @file test_colwise_reduce.cpp
 *
 * @brief Unit testing of column-wise reduction
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_reduce.h>
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

const index_t max_nrows = 32;
index_t test_nrows[] = { 1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 16, 24, 28, 32 };
const unsigned int ntest_nrows = sizeof(test_nrows) / sizeof(index_t);


#define DEFINE_COLWISE_REDUCE_CASE( Name ) \
	SIMPLE_CASE( colwise_reduce, Name ) { \
		const index_t n = 3; \
		dense_matrix<double> src(max_nrows, n); \
		fill_rand(src); \
		dense_row<double> rrow(n); \
		dense_row<double> drow1(1, zero()); \
		dense_row<double> drow(n, zero()); \
		for (unsigned k = 0; k < ntest_nrows; ++k) { \
			index_t cl = test_nrows[k]; \
			ref_col<double> s1 = src(range(0, cl), 0); \
			ref_block<double> sn = src(range(0, cl), whole()); \
			for (index_t j = 0; j < n; ++j) rrow[j] = Name(sn.column(j)); \
			colwise_##Name(s1, drow1); \
			ASSERT_APPROX( drow1[0], rrow[0], 1.0e-12 ); \
			colwise_##Name(sn, drow); \
			ASSERT_MAT_APPROX( 1, n, drow, rrow, 1.0e-12 ); } }


DEFINE_COLWISE_REDUCE_CASE( sum )
DEFINE_COLWISE_REDUCE_CASE( mean )
DEFINE_COLWISE_REDUCE_CASE( maximum )
DEFINE_COLWISE_REDUCE_CASE( minimum )

BEGIN_TPACK( colwise_reduce )
	ADD_SIMPLE_CASE( colwise_reduce, sum )
	ADD_SIMPLE_CASE( colwise_reduce, mean )
	ADD_SIMPLE_CASE( colwise_reduce, maximum )
	ADD_SIMPLE_CASE( colwise_reduce, minimum )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( colwise_reduce )
END_MAIN_SUITE

