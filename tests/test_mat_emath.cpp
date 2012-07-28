/**
 * @file test_mat_emath.cpp
 *
 * Unit test of elementary functions on matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_emath.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

const int default_m = 8;
const int default_n = 6;
const index_t LDim = 12;

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}

MN_CASE( mat_emath, pow )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n); fill_ran(A, 0.0, 5.0);
	mat_t B(m, n); fill_ran(B, 0.0, 2.0);
	double c = 1.5;
	double tol = 1.0e-12;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = std::pow(A[i], B[i]);

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = std::pow(A[i], c);

	mat_t AB = pow(A, B);
	ASSERT_TRUE( is_approx(AB, AB_r, tol) );

	mat_t AC = pow(A, c);
	ASSERT_TRUE( is_approx(AC, AC_r, tol) );
}


MN_CASE( mat_emath, floor )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::floor(A[i]);

	mat_t R = floor(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(floor(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_emath, ceil )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::ceil(A[i]);

	mat_t R = ceil(A);
	ASSERT_TRUE( is_equal(R, R_r) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(ceil(A), R_s);
	ASSERT_TRUE( is_equal(R_s, R_r) );
}


MN_CASE( mat_emath, exp )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -1.0, 3.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::exp(A[i]);

	mat_t R = exp(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(exp(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_emath, log )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, 1.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::log(A[i]);

	mat_t R = log(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(log(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_emath, log10 )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, 1.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::log10(A[i]);

	mat_t R = log10(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(log10(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_emath, sin )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::sin(A[i]);

	mat_t R = sin(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(sin(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, cos )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::cos(A[i]);

	mat_t R = cos(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(cos(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, tan )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-10;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::tan(A[i]);

	mat_t R = tan(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(tan(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_emath, asin )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -1.0, 1.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::asin(A[i]);

	mat_t R = asin(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(asin(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, acos )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -1.0, 1.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::acos(A[i]);

	mat_t R = acos(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(acos(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, atan )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::atan(A[i]);

	mat_t R = atan(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(atan(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, atan2 )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n); fill_ran(A, 1.0, 10.0);
	mat_t B(m, n); fill_ran(B, 1.0, 10.0);
	double tol = 1.0e-15;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = std::atan2(A[i], B[i]);

	mat_t AB = atan2(A, B);
	ASSERT_TRUE( is_approx(AB, AB_r, tol) );
}


MN_CASE( mat_emath, sinh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -3.0, 3.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::sinh(A[i]);

	mat_t R = sinh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(sinh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, cosh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -3.0, 3.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::cosh(A[i]);

	mat_t R = cosh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(cosh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, tanh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -5.0, 5.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::tanh(A[i]);

	mat_t R = tanh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(tanh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


#ifdef LMAT_HAS_CXX11_MATH

MN_CASE( mat_emath, round )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::round(A[i]);

	mat_t R = round(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(round(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, trunc )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::trunc(A[i]);

	mat_t R = trunc(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(trunc(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}


MN_CASE( mat_emath, cbrt )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -10.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::cbrt(A[i]);

	mat_t R = cbrt(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(cbrt(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, hypot )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	mat_t A(m, n); fill_ran(A, -5.0, 5.0);
	mat_t B(m, n); fill_ran(B, -5.0, 5.0);
	double tol = 1.0e-15;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = std::hypot(A[i], B[i]);

	mat_t AB = hypot(A, B);
	ASSERT_TRUE( is_approx(AB, AB_r, tol) );
}


MN_CASE( mat_emath, exp2 )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -1.0, 4.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::exp2(A[i]);

	mat_t R = exp2(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(exp2(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, log2 )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, 1.0, 10.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::log2(A[i]);

	mat_t R = log2(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(log2(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, expm1 )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -1.0, 1.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::expm1(A[i]);

	mat_t R = expm1(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(expm1(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, log1p )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -0.5, 1.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::log1p(A[i]);

	mat_t R = log1p(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(log1p(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, asinh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -5.0, 5.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::asinh(A[i]);

	mat_t R = asinh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(asinh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, acosh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, 1.0, 3.0);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::acosh(A[i]);

	mat_t R = acosh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(acosh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}

MN_CASE( mat_emath, atanh )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? default_m : M;
	const index_t n = N == 0 ? default_n : N;

	double tol = 1.0e-12;
	mat_t A(m, n); fill_ran(A, -0.9, 0.9);

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::atanh(A[i]);

	mat_t R = atanh(A);
	ASSERT_TRUE( is_approx(R, R_r, tol) );

	mat_t R_s(m, n, fill_value(0.0));
	linear_scalar_evaluate(atanh(A), R_s);
	ASSERT_TRUE( is_approx(R_s, R_r, tol) );
}



#endif


BEGIN_TPACK( mat_pow )
	ADD_MN_CASE_3X3( mat_emath, pow, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_floor )
	ADD_MN_CASE_3X3( mat_emath, floor, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_ceil )
	ADD_MN_CASE_3X3( mat_emath, ceil, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_exp )
	ADD_MN_CASE_3X3( mat_emath, exp, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_log )
	ADD_MN_CASE_3X3( mat_emath, log, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_log10 )
	ADD_MN_CASE_3X3( mat_emath, log10, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_sin )
	ADD_MN_CASE_3X3( mat_emath, sin, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_cos )
	ADD_MN_CASE_3X3( mat_emath, cos, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_tan )
	ADD_MN_CASE_3X3( mat_emath, tan, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_asin )
	ADD_MN_CASE_3X3( mat_emath, asin, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_acos )
	ADD_MN_CASE_3X3( mat_emath, acos, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_atan )
	ADD_MN_CASE_3X3( mat_emath, atan, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_atan2 )
	ADD_MN_CASE_3X3( mat_emath, atan2, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_sinh )
	ADD_MN_CASE_3X3( mat_emath, sinh, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_cosh )
	ADD_MN_CASE_3X3( mat_emath, cosh, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_tanh )
	ADD_MN_CASE_3X3( mat_emath, tanh, default_m, default_n )
END_TPACK

#ifdef LMAT_HAS_CXX11_MATH

BEGIN_TPACK( mat_round )
	ADD_MN_CASE_3X3( mat_emath, round, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_trunc )
	ADD_MN_CASE_3X3( mat_emath, trunc, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_cbrt )
	ADD_MN_CASE_3X3( mat_emath, cbrt, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_hypot )
	ADD_MN_CASE_3X3( mat_emath, hypot, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_exp2 )
	ADD_MN_CASE_3X3( mat_emath, exp2, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_log2 )
	ADD_MN_CASE_3X3( mat_emath, log2, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_expm1 )
	ADD_MN_CASE_3X3( mat_emath, expm1, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_log1p )
	ADD_MN_CASE_3X3( mat_emath, log1p, default_m, default_n )
END_TPACK


BEGIN_TPACK( mat_asinh )
	ADD_MN_CASE_3X3( mat_emath, asinh, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_acosh )
	ADD_MN_CASE_3X3( mat_emath, acosh, default_m, default_n )
END_TPACK

BEGIN_TPACK( mat_atanh )
	ADD_MN_CASE_3X3( mat_emath, atanh, default_m, default_n )
END_TPACK

#endif


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_floor )
	ADD_TPACK( mat_ceil )

	ADD_TPACK( mat_pow )

	ADD_TPACK( mat_exp )
	ADD_TPACK( mat_log )
	ADD_TPACK( mat_log10 )

	ADD_TPACK( mat_sin )
	ADD_TPACK( mat_cos )
	ADD_TPACK( mat_tan )

	ADD_TPACK( mat_asin )
	ADD_TPACK( mat_acos )
	ADD_TPACK( mat_atan )
	ADD_TPACK( mat_atan2 )

	ADD_TPACK( mat_sinh )
	ADD_TPACK( mat_cosh )
	ADD_TPACK( mat_tanh )

#ifdef LMAT_HAS_CXX11_MATH

	ADD_TPACK( mat_round )
	ADD_TPACK( mat_trunc )

	ADD_TPACK( mat_cbrt )
	ADD_TPACK( mat_hypot )

	ADD_TPACK( mat_exp2 )
	ADD_TPACK( mat_log2 )
	ADD_TPACK( mat_expm1 )
	ADD_TPACK( mat_log1p )

	ADD_TPACK( mat_asinh )
	ADD_TPACK( mat_acosh )
	ADD_TPACK( mat_atanh )

#endif

END_MAIN_SUITE





