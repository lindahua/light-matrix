/**
 * @file test_rowwise_reduce.cpp
 *
 * @brief Unit testing of row-wise reduction
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


#define DEFINE_ROWWISE_REDUCE_CASE( Name ) \
		SIMPLE_CASE( rowwise_reduce, Name ) { \
			const index_t n = 9; \
			dense_matrix<double> src(max_nrows, n); \
			dense_matrix<double> tsrc(n, max_nrows); \
			fill_rand(src); \
			for (index_t j = 0; j < n; ++j) { \
				for (index_t i = 0; i < max_nrows; ++i) \
					tsrc(j, i) = src(i, j); } \
			for (unsigned k = 0; k < ntest_nrows; ++k) { \
				index_t nr = test_nrows[k]; \
				auto s1 = src(range(0, nr), 1); \
				auto sn = src(range(0, nr), whole()); \
				auto t1 = tsrc(1, range(0, nr)); \
				auto tn = tsrc(whole(), range(0, nr)); \
				dense_col<double> d(nr); \
				dense_col<double> r(nr); \
				rowwise_##Name(s1, d); \
				colwise_##Name(t1, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); \
				rowwise_##Name(sn, d); \
				colwise_##Name(tn, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); } }


DEFINE_ROWWISE_REDUCE_CASE( sum )
DEFINE_ROWWISE_REDUCE_CASE( mean )
DEFINE_ROWWISE_REDUCE_CASE( maximum )
DEFINE_ROWWISE_REDUCE_CASE( minimum )

BEGIN_TPACK( rowwise_reduce )
	ADD_SIMPLE_CASE( rowwise_reduce, sum )
	ADD_SIMPLE_CASE( rowwise_reduce, mean )
	ADD_SIMPLE_CASE( rowwise_reduce, maximum )
	ADD_SIMPLE_CASE( rowwise_reduce, minimum )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( rowwise_reduce )
END_MAIN_SUITE

