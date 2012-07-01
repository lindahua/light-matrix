/**
 * @file test_mat_arith.cpp
 *
 * Unit testing for Matrix arithmetics
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/matrix_arith.h>

using namespace lmat;
using namespace lmat::test;

const int default_m = 8;
const int default_n = 6;
const index_t LDim = 12;

typedef dense_matrix<double> dmat_t;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert( is_mat_xpr<binary_ewise_expr<add_op<double>, dmat_t, dmat_t> >::value,
			"Expression interface verification failed");

	static_assert( is_mat_xpr<binary_fix1_ewise_expr<add_op<double>, dmat_t> >::value,
			"Expression interface verification failed");

	static_assert( is_mat_xpr<binary_fix2_ewise_expr<add_op<double>, dmat_t> >::value,
			"Expression interface verification failed");

	static_assert( is_linear_vector_evaluator<
			binary_ewise_linear_evaluator<add_op<double>, dmat_t, dmat_t>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_linear_vector_evaluator<
			binary_fix1_ewise_linear_evaluator<add_op<double>, dmat_t>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_linear_vector_evaluator<
			binary_fix2_ewise_linear_evaluator<add_op<double>, dmat_t>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_percol_vector_evaluator<
			binary_ewise_percol_evaluator<add_op<double>, dmat_t, dmat_t>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_percol_vector_evaluator<
			binary_fix1_ewise_percol_evaluator<add_op<double>, dmat_t>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_percol_vector_evaluator<
			binary_fix2_ewise_percol_evaluator<add_op<double>, dmat_t>, double>::value,
			"Evaluator interface verification failed");
#endif


MN_CASE( mat_arith, add )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<add_op<double>, mat_t, mat_t> AB_t;
	typedef binary_fix2_ewise_expr<add_op<double>, mat_t> AC_t;
	typedef binary_fix1_ewise_expr<add_op<double>, mat_t> CB_t;

	static_assert(is_same<decltype(A + B), AB_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(A + c), AC_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(c + B), CB_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c + B[i];

	// default evaluation

	mat_t AB = A + B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A + c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c + B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A + B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A + c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c + B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, add_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	scoped_array<double> sb(LDim * n);

	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double(i + 1);
	for (index_t i = 0; i < LDim * n; ++i) sb[i] = double(2 * i + 3);

	mat_ex_t A(sa.ptr_begin(), m, n, LDim);
	mat_ex_t B(sb.ptr_begin(), m, n, LDim);
	double c = 7.0;

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AB_r(i,j) = A(i,j) + B(i,j);

	mat_t AC_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AC_r(i,j) = A(i,j) + c;

	mat_t CB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) CB_r(i,j) = c + B(i,j);

	// default evaluation

	mat_t AB = A + B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A + c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c + B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A + B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A + c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c + B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, add_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	// default evaluation

	mat_t AB(A);
	AB += B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC += c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}



MN_CASE( mat_arith, sub )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<sub_op<double>, mat_t, mat_t> AB_t;
	typedef binary_fix2_ewise_expr<sub_op<double>, mat_t> AC_t;
	typedef binary_fix1_ewise_expr<sub_op<double>, mat_t> CB_t;

	static_assert(is_same<decltype(A - B), AB_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(A - c), AC_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(c - B), CB_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c - B[i];

	// default evaluation

	mat_t AB = A - B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A - c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c - B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A - B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A - c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c - B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, sub_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	scoped_array<double> sb(LDim * n);

	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double(i + 1);
	for (index_t i = 0; i < LDim * n; ++i) sb[i] = double(2 * i + 3);

	mat_ex_t A(sa.ptr_begin(), m, n, LDim);
	mat_ex_t B(sb.ptr_begin(), m, n, LDim);
	double c = 7.0;

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AB_r(i,j) = A(i,j) - B(i,j);

	mat_t AC_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AC_r(i,j) = A(i,j) - c;

	mat_t CB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) CB_r(i,j) = c - B(i,j);

	// default evaluation

	mat_t AB = A - B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A - c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c - B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A - B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A - c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c - B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, sub_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	// default evaluation

	mat_t AB(A);
	AB -= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC -= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, mul )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<mul_op<double>, mat_t, mat_t> AB_t;
	typedef binary_fix2_ewise_expr<mul_op<double>, mat_t> AC_t;
	typedef binary_fix1_ewise_expr<mul_op<double>, mat_t> CB_t;

	static_assert(is_same<decltype(A * B), AB_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(A * c), AC_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(c * B), CB_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c * B[i];

	// default evaluation

	mat_t AB = A * B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A * c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c * B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A * B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A * c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c * B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, mul_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	scoped_array<double> sb(LDim * n);

	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double(i + 1);
	for (index_t i = 0; i < LDim * n; ++i) sb[i] = double(2 * i + 3);

	mat_ex_t A(sa.ptr_begin(), m, n, LDim);
	mat_ex_t B(sb.ptr_begin(), m, n, LDim);
	double c = 7.0;

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AB_r(i,j) = A(i,j) * B(i,j);

	mat_t AC_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AC_r(i,j) = A(i,j) * c;

	mat_t CB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) CB_r(i,j) = c * B(i,j);

	// default evaluation

	mat_t AB = A * B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A * c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c * B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A * B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A * c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c * B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, mul_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	// default evaluation

	mat_t AB(A);
	AB *= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC *= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, div )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef binary_ewise_expr<div_op<double>, mat_t, mat_t> AB_t;
	typedef binary_fix2_ewise_expr<div_op<double>, mat_t> AC_t;
	typedef binary_fix1_ewise_expr<div_op<double>, mat_t> CB_t;

	static_assert(is_same<decltype(A / B), AB_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(A / c), AC_t>::value, "Expression type verification failed.");
	static_assert(is_same<decltype(c / B), CB_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c / B[i];

	// default evaluation

	mat_t AB = A / B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A / c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c / B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A / B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A / c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c / B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}


MN_CASE( mat_arith, div_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	scoped_array<double> sb(LDim * n);

	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double(i + 1);
	for (index_t i = 0; i < LDim * n; ++i) sb[i] = double(1 << (i % 5));

	mat_ex_t A(sa.ptr_begin(), m, n, LDim);
	mat_ex_t B(sb.ptr_begin(), m, n, LDim);
	double c = 4.0;

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AB_r(i,j) = A(i,j) / B(i,j);

	mat_t AC_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) AC_r(i,j) = A(i,j) / c;

	mat_t CB_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) CB_r(i,j) = c / B(i,j);

	// default evaluation

	mat_t AB = A / B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A / c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c / B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

	// by-scalars evaluation

	mat_t AB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A / B, AB_s);
	ASSERT_TRUE( is_equal(AB_s, AB_r) );

	mat_t AC_s(m, n, fill_value(0.0));
	evaluate_by_scalars(A / c, AC_s);
	ASSERT_TRUE( is_equal(AC_s, AC_r) );

	mat_t CB_s(m, n, fill_value(0.0));
	evaluate_by_scalars(c / B, CB_s);
	ASSERT_TRUE( is_equal(CB_s, CB_r) );
}

MN_CASE( mat_arith, div_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	// default evaluation

	mat_t AB(A);
	AB /= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC /= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}



BEGIN_TPACK( mat_arith_add )
	ADD_MN_CASE( mat_arith, add, 0, 0 )
	ADD_MN_CASE( mat_arith, add, 0, 1 )
	ADD_MN_CASE( mat_arith, add, 0, default_n )
	ADD_MN_CASE( mat_arith, add, 1, 0 )
	ADD_MN_CASE( mat_arith, add, 1, 1 )
	ADD_MN_CASE( mat_arith, add, 1, default_n )
	ADD_MN_CASE( mat_arith, add, default_m, 0 )
	ADD_MN_CASE( mat_arith, add, default_m, 1 )
	ADD_MN_CASE( mat_arith, add, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_add_ex )
	ADD_MN_CASE( mat_arith, add_ex, 0, 0 )
	ADD_MN_CASE( mat_arith, add_ex, 0, 1 )
	ADD_MN_CASE( mat_arith, add_ex, 0, default_n )
	ADD_MN_CASE( mat_arith, add_ex, 1, 0 )
	ADD_MN_CASE( mat_arith, add_ex, 1, 1 )
	ADD_MN_CASE( mat_arith, add_ex, 1, default_n )
	ADD_MN_CASE( mat_arith, add_ex, default_m, 0 )
	ADD_MN_CASE( mat_arith, add_ex, default_m, 1 )
	ADD_MN_CASE( mat_arith, add_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_add_ip )
	ADD_MN_CASE( mat_arith, add_ip, 0, 0 )
	ADD_MN_CASE( mat_arith, add_ip, 0, 1 )
	ADD_MN_CASE( mat_arith, add_ip, 0, default_n )
	ADD_MN_CASE( mat_arith, add_ip, 1, 0 )
	ADD_MN_CASE( mat_arith, add_ip, 1, 1 )
	ADD_MN_CASE( mat_arith, add_ip, 1, default_n )
	ADD_MN_CASE( mat_arith, add_ip, default_m, 0 )
	ADD_MN_CASE( mat_arith, add_ip, default_m, 1 )
	ADD_MN_CASE( mat_arith, add_ip, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_sub )
	ADD_MN_CASE( mat_arith, sub, 0, 0 )
	ADD_MN_CASE( mat_arith, sub, 0, 1 )
	ADD_MN_CASE( mat_arith, sub, 0, default_n )
	ADD_MN_CASE( mat_arith, sub, 1, 0 )
	ADD_MN_CASE( mat_arith, sub, 1, 1 )
	ADD_MN_CASE( mat_arith, sub, 1, default_n )
	ADD_MN_CASE( mat_arith, sub, default_m, 0 )
	ADD_MN_CASE( mat_arith, sub, default_m, 1 )
	ADD_MN_CASE( mat_arith, sub, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ex )
	ADD_MN_CASE( mat_arith, sub_ex, 0, 0 )
	ADD_MN_CASE( mat_arith, sub_ex, 0, 1 )
	ADD_MN_CASE( mat_arith, sub_ex, 0, default_n )
	ADD_MN_CASE( mat_arith, sub_ex, 1, 0 )
	ADD_MN_CASE( mat_arith, sub_ex, 1, 1 )
	ADD_MN_CASE( mat_arith, sub_ex, 1, default_n )
	ADD_MN_CASE( mat_arith, sub_ex, default_m, 0 )
	ADD_MN_CASE( mat_arith, sub_ex, default_m, 1 )
	ADD_MN_CASE( mat_arith, sub_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ip )
	ADD_MN_CASE( mat_arith, sub_ip, 0, 0 )
	ADD_MN_CASE( mat_arith, sub_ip, 0, 1 )
	ADD_MN_CASE( mat_arith, sub_ip, 0, default_n )
	ADD_MN_CASE( mat_arith, sub_ip, 1, 0 )
	ADD_MN_CASE( mat_arith, sub_ip, 1, 1 )
	ADD_MN_CASE( mat_arith, sub_ip, 1, default_n )
	ADD_MN_CASE( mat_arith, sub_ip, default_m, 0 )
	ADD_MN_CASE( mat_arith, sub_ip, default_m, 1 )
	ADD_MN_CASE( mat_arith, sub_ip, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_mul )
	ADD_MN_CASE( mat_arith, mul, 0, 0 )
	ADD_MN_CASE( mat_arith, mul, 0, 1 )
	ADD_MN_CASE( mat_arith, mul, 0, default_n )
	ADD_MN_CASE( mat_arith, mul, 1, 0 )
	ADD_MN_CASE( mat_arith, mul, 1, 1 )
	ADD_MN_CASE( mat_arith, mul, 1, default_n )
	ADD_MN_CASE( mat_arith, mul, default_m, 0 )
	ADD_MN_CASE( mat_arith, mul, default_m, 1 )
	ADD_MN_CASE( mat_arith, mul, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ex )
	ADD_MN_CASE( mat_arith, mul_ex, 0, 0 )
	ADD_MN_CASE( mat_arith, mul_ex, 0, 1 )
	ADD_MN_CASE( mat_arith, mul_ex, 0, default_n )
	ADD_MN_CASE( mat_arith, mul_ex, 1, 0 )
	ADD_MN_CASE( mat_arith, mul_ex, 1, 1 )
	ADD_MN_CASE( mat_arith, mul_ex, 1, default_n )
	ADD_MN_CASE( mat_arith, mul_ex, default_m, 0 )
	ADD_MN_CASE( mat_arith, mul_ex, default_m, 1 )
	ADD_MN_CASE( mat_arith, mul_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ip )
	ADD_MN_CASE( mat_arith, mul_ip, 0, 0 )
	ADD_MN_CASE( mat_arith, mul_ip, 0, 1 )
	ADD_MN_CASE( mat_arith, mul_ip, 0, default_n )
	ADD_MN_CASE( mat_arith, mul_ip, 1, 0 )
	ADD_MN_CASE( mat_arith, mul_ip, 1, 1 )
	ADD_MN_CASE( mat_arith, mul_ip, 1, default_n )
	ADD_MN_CASE( mat_arith, mul_ip, default_m, 0 )
	ADD_MN_CASE( mat_arith, mul_ip, default_m, 1 )
	ADD_MN_CASE( mat_arith, mul_ip, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_div )
	ADD_MN_CASE( mat_arith, div, 0, 0 )
	ADD_MN_CASE( mat_arith, div, 0, 1 )
	ADD_MN_CASE( mat_arith, div, 0, default_n )
	ADD_MN_CASE( mat_arith, div, 1, 0 )
	ADD_MN_CASE( mat_arith, div, 1, 1 )
	ADD_MN_CASE( mat_arith, div, 1, default_n )
	ADD_MN_CASE( mat_arith, div, default_m, 0 )
	ADD_MN_CASE( mat_arith, div, default_m, 1 )
	ADD_MN_CASE( mat_arith, div, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_div_ex )
	ADD_MN_CASE( mat_arith, div_ex, 0, 0 )
	ADD_MN_CASE( mat_arith, div_ex, 0, 1 )
	ADD_MN_CASE( mat_arith, div_ex, 0, default_n )
	ADD_MN_CASE( mat_arith, div_ex, 1, 0 )
	ADD_MN_CASE( mat_arith, div_ex, 1, 1 )
	ADD_MN_CASE( mat_arith, div_ex, 1, default_n )
	ADD_MN_CASE( mat_arith, div_ex, default_m, 0 )
	ADD_MN_CASE( mat_arith, div_ex, default_m, 1 )
	ADD_MN_CASE( mat_arith, div_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_div_ip )
	ADD_MN_CASE( mat_arith, div_ip, 0, 0 )
	ADD_MN_CASE( mat_arith, div_ip, 0, 1 )
	ADD_MN_CASE( mat_arith, div_ip, 0, default_n )
	ADD_MN_CASE( mat_arith, div_ip, 1, 0 )
	ADD_MN_CASE( mat_arith, div_ip, 1, 1 )
	ADD_MN_CASE( mat_arith, div_ip, 1, default_n )
	ADD_MN_CASE( mat_arith, div_ip, default_m, 0 )
	ADD_MN_CASE( mat_arith, div_ip, default_m, 1 )
	ADD_MN_CASE( mat_arith, div_ip, default_m, default_n )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_arith_add )
	ADD_TPACK( mat_arith_add_ex )
	ADD_TPACK( mat_arith_add_ip )

	ADD_TPACK( mat_arith_sub )
	ADD_TPACK( mat_arith_sub_ex )
	ADD_TPACK( mat_arith_sub_ip )

	ADD_TPACK( mat_arith_mul )
	ADD_TPACK( mat_arith_mul_ex )
	ADD_TPACK( mat_arith_mul_ip )

	ADD_TPACK( mat_arith_div )
	ADD_TPACK( mat_arith_div_ex )
	ADD_TPACK( mat_arith_div_ip )
END_MAIN_SUITE



