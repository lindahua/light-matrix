/**
 * @file test_dense_cast.cpp
 *
 * Unit testing of casting for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

using namespace lmat;
using namespace lmat::test;

#ifdef LMAT_USE_STATIC_ASSERT
static_assert( is_implicitly_convertible<float, double>::value == true,
		"implicitly_convertible test failed");
static_assert( is_implicitly_convertible<double, float>::value == false,
		"implicitly_convertible test failed");
#endif


MN_CASE( dense_cast, implicit_cast )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<float,  M, N> A(m, n);
	dense_matrix<double, M, N> Br(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		index_t v = i + 2;
		A[i] = float(v);
		Br[i] = double(v);
	}

	dense_matrix<double, M, N> B = A;

	ASSERT_TRUE( is_equal(B, Br) );
}

MN_CASE( dense_cast, explicit_cast )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> A(m, n);
	dense_matrix<float,  M, N> Br(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		int v = i + 2;
		A[i] = double(v);
		Br[i] = float(v);
	}

	dense_matrix<float, M, N> B1 = A.cast<float>();
	dense_matrix<float, M, N> B2 = cast(A, type<float>());

	ASSERT_TRUE( is_equal(B1, Br) );
	ASSERT_TRUE( is_equal(B2, Br) );
}


MN_CASE( ref_cast, implicit_cast )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<float> sA(m * n);
	scoped_array<double> sB(m * n);

	ref_matrix<float,  M, N> A(sA.ptr_begin(), m, n);
	dense_matrix<double, M, N> Br(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		index_t v = i + 2;
		A[i] = float(v);
		Br[i] = double(v);
	}

	ref_matrix<double, M, N> B(sB.ptr_begin(), m, n);
	B = A;

	ASSERT_TRUE( is_equal(B, Br) );
}

MN_CASE( ref_cast, explicit_cast )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<double> sA(m * n);
	scoped_array<float> sB(m * n);

	ref_matrix<double, M, N> A(sA.ptr_begin(), m, n);
	dense_matrix<float, M, N> Br(m, n);

	for (index_t i = 0; i < m * n; ++i)
	{
		index_t v = i + 2;
		A[i] = double(v);
		Br[i] = float(v);
	}

	ref_matrix<float, M, N> B(sB.ptr_begin(), m, n);
	B = A.cast<float>();

	ASSERT_TRUE( is_equal(B, Br) );
}


MN_CASE( ref_ex_cast, implicit_cast )
{
	const index_t ldim_a = 7;
	const index_t ldim_b = 8;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<float> sA(ldim_a * n);
	scoped_array<double> sB(ldim_b * n);

	ref_matrix_ex<float,  M, N> A(sA.ptr_begin(), m, n, ldim_a);
	dense_matrix<double, M, N> Br(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			index_t v = i + j + 2;
			A(i, j) = float(v);
			Br(i, j) = double(v);
		}
	}

	ref_matrix_ex<double, M, N> B(sB.ptr_begin(), m, n, ldim_b);
	B = A;

	ASSERT_TRUE( is_equal(B, Br) );
}

MN_CASE( ref_ex_cast, explicit_cast )
{
	const index_t ldim_a = 7;
	const index_t ldim_b = 8;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<double> sA(ldim_a * n);
	scoped_array<float> sB(ldim_b * n);

	ref_matrix_ex<double,  M, N> A(sA.ptr_begin(), m, n, ldim_a);
	dense_matrix<float, M, N> Br(m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			index_t v = i + j + 2;
			A(i, j) = double(v);
			Br(i, j) = float(v);
		}
	}

	ref_matrix_ex<float, M, N> B(sB.ptr_begin(), m, n, ldim_b);
	B = A.cast<float>();

	ASSERT_TRUE( is_equal(B, Br) );
}



BEGIN_TPACK( dense_implicit_cast )
	ADD_MN_CASE_3X3( dense_cast, implicit_cast, 5, 6 );
END_TPACK

BEGIN_TPACK( dense_explicit_cast )
	ADD_MN_CASE_3X3( dense_cast, explicit_cast, 5, 6 );
END_TPACK

BEGIN_TPACK( ref_implicit_cast )
	ADD_MN_CASE_3X3( ref_cast, implicit_cast, 5, 6 );
END_TPACK

BEGIN_TPACK( ref_explicit_cast )
	ADD_MN_CASE_3X3( ref_cast, explicit_cast, 5, 6 );
END_TPACK

BEGIN_TPACK( ref_ex_implicit_cast )
	ADD_MN_CASE_3X3( ref_ex_cast, implicit_cast, 5, 6 );
END_TPACK

BEGIN_TPACK( ref_ex_explicit_cast )
	ADD_MN_CASE_3X3( ref_ex_cast, explicit_cast, 5, 6 );
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( dense_implicit_cast )
	ADD_TPACK( dense_explicit_cast )
	ADD_TPACK( ref_implicit_cast )
	ADD_TPACK( ref_explicit_cast )
	ADD_TPACK( ref_ex_implicit_cast )
	ADD_TPACK( ref_ex_explicit_cast )
END_MAIN_SUITE





