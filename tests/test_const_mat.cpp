/**
 * @file test_const_mat.cpp
 *
 * Unit testing of const_matrix
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/const_matrix.h>
#include <light_mat/matrix/dense_matrix.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::const_matrix<double, 0, 0>;
template class lmat::const_matrix<double, 0, 4>;
template class lmat::const_matrix<double, 3, 0>;
template class lmat::const_matrix<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT
static_assert(lmat::is_mat_xpr<lmat::const_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::const_matrix<double> >::value, "Interface verification failed.");
#endif

MN_CASE( const_mat, constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const double val = 12.5;

	const_matrix<double, M, N> a(m, n, val);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.size(), size_t(m * n));

	ASSERT_EQ(a.value(), val);
}

MN_CASE( const_mat, access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const double val = 12.5;

	const_matrix<double, M, N> a(m, n, val);
	dense_matrix<double, M, N> r(m, n, fill(val));

	ASSERT_MAT_EQ(m, n, a, r);
}


MN_CASE( const_mat, evaluates )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const double val = 12.5;

	const_matrix<double, M, N> a(m, n, val);

	dense_matrix<double, M, N> b = a;
	dense_matrix<double, M, N> r(m, n, fill(val));

	ASSERT_MAT_EQ(m, n, b, r);
}


BEGIN_TPACK( const_mat_constructs )
	ADD_MN_CASE_3X3( const_mat, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( const_mat_access )
	ADD_MN_CASE_3X3( const_mat, access, 3, 4 )
END_TPACK

BEGIN_TPACK( const_mat_evaluates )
	ADD_MN_CASE_3X3( const_mat, evaluates, 3, 4 )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( const_mat_constructs )
	ADD_TPACK( const_mat_access )
	ADD_TPACK( const_mat_evaluates )
END_MAIN_SUITE





