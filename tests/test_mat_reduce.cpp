/**
 * @file test_mat_reduce.cpp
 *
 * Unit testing for matrix reduction
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/matrix_reduce.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 6;
const index_t DN = 7;
const index_t LDim = 8;

MN_CASE( mat_reduce, sum )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += A[i];

	// test

	ASSERT_EQ( reduce_by_scalars(sum_fun<double>(), A), r0 );
	ASSERT_EQ( sum(A), r0 );
}

MN_CASE( mat_reduce, sum_ex )
{
	typedef ref_matrix_ex<double, M, N> mat_ex;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	scoped_array<double> src(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) src[i] = double(i + 1);

	mat_ex A(src.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) r0 += A(i, j);

	// test

	ASSERT_EQ( reduce_by_scalars(sum_fun<double>(), A), r0 );
	ASSERT_EQ( sum(A), r0 );
}


MN_CASE( mat_reduce, mean )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += A[i];
	r0 /= double(m * n);

	// test

	ASSERT_TRUE( is_approx_scalar(mean(A), r0, 1.0e-12 ) );
}


MN_CASE( mat_reduce, maximum )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (10 - 2 * i));

	// prepare ground-truth

	double r0 = A[0];
	for (index_t i = 0; i < m * n; ++i)
	{
		if (A[i] > r0) r0 = A[i];
	}

	// test

	ASSERT_EQ( reduce_by_scalars(maximum_fun<double>(), A), r0 );
	ASSERT_EQ( maximum(A), r0 );
}

MN_CASE( mat_reduce, maximum_ex )
{
	typedef ref_matrix_ex<double, M, N> mat_ex;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	scoped_array<double> src(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) src[i] = double((i + 1) * (10 - 2 * i));

	mat_ex A(src.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	double r0 = A(0, 0);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i)
			if (A(i, j) > r0) r0 = A(i, j);

	// test

	ASSERT_EQ( reduce_by_scalars(maximum_fun<double>(), A), r0 );
	ASSERT_EQ( maximum(A), r0 );
}


MN_CASE( mat_reduce, minimum )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (i - 10));

	// prepare ground-truth

	double r0 = A[0];
	for (index_t i = 0; i < m * n; ++i)
	{
		if (A[i] < r0) r0 = A[i];
	}

	// test

	ASSERT_EQ( reduce_by_scalars(minimum_fun<double>(), A), r0 );
	ASSERT_EQ( minimum(A), r0 );
}

MN_CASE( mat_reduce, minimum_ex )
{
	typedef ref_matrix_ex<double, M, N> mat_ex;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	scoped_array<double> src(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) src[i] = double((i + 1) * (i - 10));

	mat_ex A(src.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	double r0 = A(0, 0);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i)
			if (A(i, j) < r0) r0 = A(i, j);

	// test

	ASSERT_EQ( reduce_by_scalars(minimum_fun<double>(), A), r0 );
	ASSERT_EQ( minimum(A), r0 );
}


MN_CASE( mat_reduce, dot )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(i + 2);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += A[i] * B[i];

	// test

	ASSERT_EQ( dot(A, B), r0 );
}


MN_CASE( mat_reduce, L1norm )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i - 10);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += std::abs(A[i]);

	// test

	ASSERT_EQ( L1norm(A), r0 );
}

MN_CASE( mat_reduce, sqL2norm )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i - 10);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += A[i] * A[i];

	// test

	ASSERT_EQ( sqL2norm(A), r0 );
}


MN_CASE( mat_reduce, L2norm )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i - 10);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i) r0 += A[i] * A[i];
	r0 = std::sqrt(r0);

	// test

	ASSERT_EQ( L2norm(A), r0 );
}


MN_CASE( mat_reduce, Linfnorm )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i - 10);

	// prepare ground-truth

	double r0 = 0.0;
	for (index_t i = 0; i < m * n; ++i)
	{
		if (std::abs(A[i]) > r0) r0 = std::abs(A[i]);
	}

	// test

	ASSERT_EQ( Linfnorm(A), r0 );
}


MN_CASE( mat_reduce, nrmdot )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(i + 2);

	// prepare ground-truth

	double raa = 0.0;
	double rab = 0.0;
	double rbb = 0.0;

	for (index_t i = 0; i < m * n; ++i)
	{
		rab += A[i] * B[i];
		raa += A[i] * A[i];
		rbb += B[i] * B[i];
	}

	double r0 = rab / (std::sqrt(raa) * std::sqrt(rbb));
	ASSERT_TRUE( r0 <= 1.0 && r0 >= -1.0 );

	// test

	ASSERT_EQ( nrmdot(A, B), r0 );
}


BEGIN_TPACK( mat_sum )
	ADD_MN_CASE_3X3( mat_reduce, sum, DM, DN )
END_TPACK

BEGIN_TPACK( mat_sum_ex )
	ADD_MN_CASE_3X3( mat_reduce, sum, DM, DN )
END_TPACK

BEGIN_TPACK( mat_mean )
	ADD_MN_CASE_3X3( mat_reduce, mean, DM, DN )
END_TPACK

BEGIN_TPACK( mat_maximum )
	ADD_MN_CASE_3X3( mat_reduce, maximum, DM, DN )
END_TPACK

BEGIN_TPACK( mat_maximum_ex )
	ADD_MN_CASE_3X3( mat_reduce, maximum_ex, DM, DN )
END_TPACK

BEGIN_TPACK( mat_minimum )
	ADD_MN_CASE_3X3( mat_reduce, minimum, DM, DN )
END_TPACK

BEGIN_TPACK( mat_minimum_ex )
	ADD_MN_CASE_3X3( mat_reduce, minimum_ex, DM, DN )
END_TPACK

BEGIN_TPACK( mat_dot )
	ADD_MN_CASE_3X3( mat_reduce, dot, DM, DN )
END_TPACK

BEGIN_TPACK( mat_L1norm )
	ADD_MN_CASE_3X3( mat_reduce, L1norm, DM, DN )
END_TPACK

BEGIN_TPACK( mat_sqL2norm )
	ADD_MN_CASE_3X3( mat_reduce, sqL2norm, DM, DN )
END_TPACK

BEGIN_TPACK( mat_L2norm )
	ADD_MN_CASE_3X3( mat_reduce, L2norm, DM, DN )
END_TPACK

BEGIN_TPACK( mat_Linfnorm )
	ADD_MN_CASE_3X3( mat_reduce, Linfnorm, DM, DN )
END_TPACK

BEGIN_TPACK( mat_nrmdot )
	ADD_MN_CASE_3X3( mat_reduce, nrmdot, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_sum )
	ADD_TPACK( mat_sum_ex )
	ADD_TPACK( mat_mean )
	ADD_TPACK( mat_maximum )
	ADD_TPACK( mat_maximum_ex )
	ADD_TPACK( mat_minimum )
	ADD_TPACK( mat_minimum_ex )

	ADD_TPACK( mat_dot )
	ADD_TPACK( mat_L1norm )
	ADD_TPACK( mat_sqL2norm )
	ADD_TPACK( mat_L2norm )
	ADD_TPACK( mat_Linfnorm )
	ADD_TPACK( mat_nrmdot )
END_MAIN_SUITE




