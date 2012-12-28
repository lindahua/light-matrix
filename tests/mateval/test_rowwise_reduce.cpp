/**
 * @file test_rowwise_reduce.cpp
 *
 * @brief Unit testing of row-wise reduction
 *
 * @author Dahua Lin
 */


#include "../test_base.h"
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
				auto s1 = src(range(0, nr), 0); \
				auto sn = src(range(0, nr), whole()); \
				auto t1 = tsrc(range(0, 1), range(0, nr)); \
				auto tn = tsrc(whole(), range(0, nr)); \
				dense_col<double> d(nr); \
				dense_col<double> r(nr); \
				rowwise_##Name(s1, d); \
				colwise_##Name(t1, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); \
				rowwise_##Name(sn, d); \
				colwise_##Name(tn, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); } }

#define DEFINE_ROWWISE_REDUCE_CASE_2( Name ) \
		SIMPLE_CASE( rowwise_reduce, Name ) { \
			const index_t n = 9; \
			dense_matrix<double> src1(max_nrows, n); \
			dense_matrix<double> src2(max_nrows, n); \
			dense_matrix<double> tsrc1(n, max_nrows); \
			dense_matrix<double> tsrc2(n, max_nrows); \
			fill_rand(src1); \
			fill_rand(src2); \
			for (index_t j = 0; j < n; ++j) { \
				for (index_t i = 0; i < max_nrows; ++i) { \
					tsrc1(j, i) = src1(i, j); \
					tsrc2(j, i) = src2(i, j); } } \
			for (unsigned k = 0; k < ntest_nrows; ++k) { \
				index_t nr = test_nrows[k]; \
				auto s11 = src1(range(0, nr), 0); \
				auto s12 = src2(range(0, nr), 0); \
				auto sn1 = src1(range(0, nr), whole()); \
				auto sn2 = src2(range(0, nr), whole()); \
				auto t11 = tsrc1(range(0, 1), range(0, nr)); \
				auto t12 = tsrc2(range(0, 1), range(0, nr)); \
				auto tn1 = tsrc1(whole(), range(0, nr)); \
				auto tn2 = tsrc2(whole(), range(0, nr)); \
				dense_col<double> d(nr); \
				dense_col<double> r(nr); \
				rowwise_##Name(s11, s12, d); \
				colwise_##Name(t11, t12, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); \
				rowwise_##Name(sn1, sn2, d); \
				colwise_##Name(tn1, tn2, r); \
				ASSERT_MAT_APPROX(nr, 1, r, d, 1.0e-12); } }


DEFINE_ROWWISE_REDUCE_CASE( sum )
DEFINE_ROWWISE_REDUCE_CASE( mean )
DEFINE_ROWWISE_REDUCE_CASE( maximum )
DEFINE_ROWWISE_REDUCE_CASE( minimum )

DEFINE_ROWWISE_REDUCE_CASE( asum )
DEFINE_ROWWISE_REDUCE_CASE( amean )
DEFINE_ROWWISE_REDUCE_CASE( amax )
DEFINE_ROWWISE_REDUCE_CASE( sqsum )

DEFINE_ROWWISE_REDUCE_CASE_2( diff_asum )
DEFINE_ROWWISE_REDUCE_CASE_2( diff_amean )
DEFINE_ROWWISE_REDUCE_CASE_2( diff_amax )
DEFINE_ROWWISE_REDUCE_CASE_2( diff_sqsum )

DEFINE_ROWWISE_REDUCE_CASE_2( dot )

BEGIN_TPACK( rowwise_reduce )
	ADD_SIMPLE_CASE( rowwise_reduce, sum )
	ADD_SIMPLE_CASE( rowwise_reduce, mean )
	ADD_SIMPLE_CASE( rowwise_reduce, maximum )
	ADD_SIMPLE_CASE( rowwise_reduce, minimum )

	ADD_SIMPLE_CASE( rowwise_reduce, asum )
	ADD_SIMPLE_CASE( rowwise_reduce, amean )
	ADD_SIMPLE_CASE( rowwise_reduce, amax )
	ADD_SIMPLE_CASE( rowwise_reduce, sqsum )

	ADD_SIMPLE_CASE( rowwise_reduce, diff_asum )
	ADD_SIMPLE_CASE( rowwise_reduce, diff_amean )
	ADD_SIMPLE_CASE( rowwise_reduce, diff_amax )
	ADD_SIMPLE_CASE( rowwise_reduce, diff_sqsum )

	ADD_SIMPLE_CASE( rowwise_reduce, dot )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( rowwise_reduce )
END_MAIN_SUITE

