/**
 * @file test_mat_arith.cpp
 *
 * Unit testing for Basic matrix arithmetics
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
	static_assert( is_mat_xpr<unary_ewise_expr<neg_op<double>, dmat_t> >::value,
		"Expression interface verification failed");

	static_assert( is_mat_xpr<binary_ewise_expr<add_op<double>, dmat_t, dmat_t> >::value,
			"Expression interface verification failed");

	static_assert( is_linear_vector_evaluator<
			binary_ewise_linear_evaluator<add_op<double>, dmat_t, dmat_t, false, false>, double>::value,
			"Evaluator interface verification failed");

	static_assert( is_percol_vector_evaluator<
			binary_ewise_percol_evaluator<add_op<double>, dmat_t, dmat_t, false, false>, double>::value,
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


MN_CASE( mat_arith, neg )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<neg_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(-A), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = -A[i];

	// default evaluation

	mat_t R = -A;
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(-A, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, neg_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double((i+1) * (i % 3 - 1));
	mat_ex_t A(sa.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) R_r(i,j) = -A(i,j);

	// default evaluation

	mat_t R = -A;
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(-A, R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, abs )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<abs_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(abs(A)), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::abs(A[i]);

	// default evaluation

	mat_t R = abs(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(abs(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, abs_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double((i+1) * (i % 3 - 1));
	mat_ex_t A(sa.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) R_r(i,j) = std::abs(A(i,j));

	// default evaluation

	mat_t R = abs(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(abs(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, sqr )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<sqr_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(sqr(A)), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = A[i] * A[i];

	// default evaluation

	mat_t R = sqr(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(sqr(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, sqr_ex )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef ref_matrix_ex<double, M, N> mat_ex_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	scoped_array<double> sa(LDim * n);
	for (index_t i = 0; i < LDim * n; ++i) sa[i] = double((i+1) * (i % 3 - 1));
	mat_ex_t A(sa.ptr_begin(), m, n, LDim);

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t j = 0; j < n; ++j)
		for (index_t i = 0; i < m; ++i) R_r(i,j) = A(i,j) * A(i,j);

	// default evaluation

	mat_t R = sqr(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(sqr(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_arith, sqrt )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;
	const double tol = 1.0e-10;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i+2);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<sqrt_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(sqrt(A)), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::sqrt(A[i]);

	// default evaluation

	mat_t R = sqrt(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(sqrt(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_arith, rcp )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;
	const double tol = 1.0e-10;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i+2);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<rcp_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(rcp(A)), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = 1.0 / A[i];

	// default evaluation

	mat_t R = rcp(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(rcp(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_arith, rsqrt )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;
	const double tol = 1.0e-10;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i+2);

	// type verification

#ifdef LMAT_HAS_DECLTYPE
	typedef unary_ewise_expr<rsqrt_op<double>, mat_t> R_t;
	static_assert(is_same<decltype(rsqrt(A)), R_t>::value, "Expression type verification failed.");
#endif

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = 1.0 / std::sqrt(A[i]);

	// default evaluation

	mat_t R = rsqrt(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	// by-scalars evaluation

	mat_t R_s(m, n, fill_value(0.0));
	evaluate_by_scalars(rsqrt(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}



BEGIN_TPACK( mat_arith_add )
	ADD_MN_CASE_3X3( mat_arith, add, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_add_ex )
	ADD_MN_CASE_3X3( mat_arith, add_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_add_ip )
	ADD_MN_CASE_3X3( mat_arith, add_ip, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_sub )
	ADD_MN_CASE_3X3( mat_arith, sub, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ex )
	ADD_MN_CASE_3X3( mat_arith, sub_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ip )
	ADD_MN_CASE_3X3( mat_arith, sub_ip, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_mul )
	ADD_MN_CASE_3X3( mat_arith, mul, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ex )
	ADD_MN_CASE_3X3( mat_arith, mul_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ip )
	ADD_MN_CASE_3X3( mat_arith, mul_ip, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_div )
	ADD_MN_CASE_3X3( mat_arith, div, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_div_ex )
	ADD_MN_CASE_3X3( mat_arith, div_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_div_ip )
	ADD_MN_CASE_3X3( mat_arith, div_ip, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_neg )
	ADD_MN_CASE_3X3( mat_arith, neg, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_neg_ex )
	ADD_MN_CASE_3X3( mat_arith, neg_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_abs )
	ADD_MN_CASE_3X3( mat_arith, abs, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_arith_abs_ex )
	ADD_MN_CASE_3X3( mat_arith, abs_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sqr )
	ADD_MN_CASE_3X3( mat_arith, sqr, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sqr_ex )
	ADD_MN_CASE_3X3( mat_arith, sqr_ex, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_sqrt )
	ADD_MN_CASE_3X3( mat_arith, sqrt, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_rcp )
	ADD_MN_CASE_3X3( mat_arith, rcp, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_arith_rsqrt )
	ADD_MN_CASE_3X3( mat_arith, rsqrt, default_m, default_n )
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

	ADD_TPACK( mat_arith_neg )
	ADD_TPACK( mat_arith_neg_ex )

	ADD_TPACK( mat_arith_abs )
	ADD_TPACK( mat_arith_abs_ex )

	ADD_TPACK( mat_arith_sqr )
	ADD_TPACK( mat_arith_sqr_ex )

	ADD_TPACK( mat_arith_sqrt )

	ADD_TPACK( mat_arith_rcp )

	ADD_TPACK( mat_arith_rsqrt )
END_MAIN_SUITE



