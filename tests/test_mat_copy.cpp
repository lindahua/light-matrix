/**
 * @file test_mat_copy.cpp
 *
 * Unit testing of matrix copy
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_matrix_ex.h>
#include <light_mat/core/array.h>

using namespace lmat;
using namespace lmat::test;

MN_CASE( mat_copy, copy )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	darray<double> src(m * n);
	darray<double> dst(m * n);
	fill(dst, 0.0);

	for (index_t i = 0; i < m * n; ++i) src[i] = double(i + 2);

	ref_matrix<double, M, N> a(src.ptr_begin(), m, n);
	ref_matrix<double, M, N> b(dst.ptr_begin(), m, n);

	copy(a, b);

	ASSERT_MAT_EQ(m, n, a, b);
	ASSERT_VEC_EQ(m * n, src, dst);
}


MN_CASE( mat_copy, copy_ex )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	const index_t ldim_a = 7;
	const index_t ldim_b = 8;

	darray<double> src(ldim_a * n);
	darray<double> dst(ldim_b * n);
	fill(dst, 0.0);

	for (index_t i = 0; i < ldim_a * n; ++i) src[i] = double(i + 2);

	ref_matrix_ex<double, M, N> a(src.ptr_begin(), m, n, ldim_a);
	ref_matrix_ex<double, M, N> b(dst.ptr_begin(), m, n, ldim_b);

	copy(a, b);

	ASSERT_MAT_EQ(m, n, a, b);

	darray<double> ref(ldim_b * n);
	fill(ref, 0.0);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim_b] = src[i + j * ldim_a];
	}

	ASSERT_VEC_EQ(ldim_b * n, ref, dst);
}


BEGIN_TPACK( mat_copy )
	ADD_MN_CASE_3X3( mat_copy, copy, 5, 6 );
END_TPACK

BEGIN_TPACK( mat_copy_ex )
	ADD_MN_CASE_3X3( mat_copy, copy_ex, 5, 6 );
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_copy )
	ADD_TPACK( mat_copy_ex )
END_MAIN_SUITE









