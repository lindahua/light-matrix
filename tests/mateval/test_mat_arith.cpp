/**
 * @file test_mat_arith.cpp
 *
 * Unit testing for Basic matrix arithmetics
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/mateval/mat_arith.h>

using namespace lmat;
using namespace lmat::test;

typedef dense_matrix<double> dmat_t;


template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}


template<class A, class B>
bool my_is_equal(const A& a, const B& b)
{
	if ( have_same_shape(a, b) )
	{
		index_t m = a.nrows();
		index_t n = a.ncolumns();

		return ltest::test_matrix_equal(m, n, a, b);
	}
	else
	{
		return false;
	}
}

MN_CASE( mat_arith, add )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = A + c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = c + B;
	ASSERT_TRUE( my_is_equal(CB, CB_r) );

}


MN_CASE( mat_arith, add_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC(A);
	AC += c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );
}



MN_CASE( mat_arith, sub )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = A - c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = c - B;
	ASSERT_TRUE( my_is_equal(CB, CB_r) );

}

MN_CASE( mat_arith, sub_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC(A);
	AC -= c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, mul )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = A * c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = c * B;
	ASSERT_TRUE( my_is_equal(CB, CB_r) );
}

MN_CASE( mat_arith, mul_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC(A);
	AC *= c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, div )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = A / c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = c / B;
	ASSERT_TRUE( my_is_equal(CB, CB_r) );
}


MN_CASE( mat_arith, div_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

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
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC(A);
	AC /= c;
	ASSERT_TRUE( my_is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, neg )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = -A[i];

	// default evaluation

	mat_t R = -A;
	ASSERT_TRUE( my_is_equal(R, R_r) );
}


MN_CASE( mat_arith, abs )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::abs(A[i]);

	// default evaluation

	mat_t R = abs(A);
	ASSERT_TRUE( my_is_equal(R, R_r) );
}


MN_CASE( mat_arith, sqr )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1));

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = A[i] * A[i];

	mat_t R = sqr(A);
	ASSERT_TRUE( my_is_equal(R, R_r) );
}


MN_CASE( mat_arith, cube )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1));

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = A[i] * A[i] * A[i];

	mat_t R = cube(A);
	ASSERT_TRUE( my_is_equal(R, R_r) );
}


MN_CASE( mat_arith, max )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n); fill_ran(A, 0.0, 10.0);
	mat_t B(m, n); fill_ran(B, 0.0, 10.0);
	double c = 5.0;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] > B[i] ? A[i] : B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] > c ? A[i] : c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c > B[i] ? c : B[i];

	mat_t AB = (max)(A, B);
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = (max)(A, c);
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = (max)(c, B);
	ASSERT_TRUE( my_is_equal(CB, CB_r) );
}

MN_CASE( mat_arith, min )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n); fill_ran(A, 0.0, 10.0);
	mat_t B(m, n); fill_ran(B, 0.0, 10.0);
	double c = 5.0;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] < B[i] ? A[i] : B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] < c ? A[i] : c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c < B[i] ? c : B[i];

	mat_t AB = (min)(A, B);
	ASSERT_TRUE( my_is_equal(AB, AB_r) );

	mat_t AC = (min)(A, c);
	ASSERT_TRUE( my_is_equal(AC, AC_r) );

	mat_t CB = (min)(c, B);
	ASSERT_TRUE( my_is_equal(CB, CB_r) );
}



BEGIN_TPACK( mat_arith_add )
	ADD_MN_CASE_3X3( mat_arith, add, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_add_ip )
	ADD_MN_CASE_3X3( mat_arith, add_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_sub )
	ADD_MN_CASE_3X3( mat_arith, sub, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ip )
	ADD_MN_CASE_3X3( mat_arith, sub_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_mul )
	ADD_MN_CASE_3X3( mat_arith, mul, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ip )
	ADD_MN_CASE_3X3( mat_arith, mul_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_div )
	ADD_MN_CASE_3X3( mat_arith, div, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_div_ip )
	ADD_MN_CASE_3X3( mat_arith, div_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_neg )
	ADD_MN_CASE_3X3( mat_arith, neg, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_abs )
	ADD_MN_CASE_3X3( mat_arith, abs, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_sqr )
	ADD_MN_CASE_3X3( mat_arith, sqr, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_cube )
	ADD_MN_CASE_3X3( mat_arith, cube, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_max )
	ADD_MN_CASE_3X3( mat_arith, max, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_min )
	ADD_MN_CASE_3X3( mat_arith, min, DM, DN )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( mat_arith_add )
	ADD_TPACK( mat_arith_add_ip )
	ADD_TPACK( mat_arith_sub )
	ADD_TPACK( mat_arith_sub_ip )
	ADD_TPACK( mat_arith_mul )
	ADD_TPACK( mat_arith_mul_ip )
	ADD_TPACK( mat_arith_div )
	ADD_TPACK( mat_arith_div_ip )

	ADD_TPACK( mat_arith_neg )
	ADD_TPACK( mat_arith_abs )
	ADD_TPACK( mat_arith_sqr )
	ADD_TPACK( mat_arith_cube )

	ADD_TPACK( mat_arith_max )
	ADD_TPACK( mat_arith_min )
END_MAIN_SUITE



