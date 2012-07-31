/*
 * @file test_gemv.cpp
 *
 * Test of generic matrix-vector multiplication
 *
 * @author Dahua Lin
 */


#include "linalg_test_base.h"

#include <light_mat/matrix/matrix_arith.h>
#include <light_mat/matrix/matrix_ewise_eval.h>

#include <light_mat/linalg/matrix_blasl2.h>

using namespace lmat;
using namespace lmat::test;


template<typename T> struct tolvalue;

template<> struct tolvalue<float>
{
	static float get() { return 2.0e-5f; }
};


template<> struct tolvalue<double>
{
	static double get() { return 1.0e-12; }
};


const int DM = 6;
const int DN = 7;
const index_t LDim = 8;


TMN_CASE( gemv_n, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;
	const T alpha = T(1.5);
	const T beta = T(1.2);
	T tol = tolvalue<T>::get();

	dense_matrix<T, M, N> A(m, n);
	dense_matrix<T, N, 1> x(n, 1);
	dense_matrix<T, M, 1> y(m, 1, zeros<T>());
	dense_matrix<T, M, 1> r(m, 1, zeros<T>());
	dense_matrix<T, M, 1> t(m, 1, zeros<T>());

	fill_rand(A, T(-1), T(1));
	fill_rand(x, T(-1), T(1));

	// y = A * x

	gemv_n(A, x, y);
	naive_mtimes(A, x, t);
	ASSERT_MAT_APPROX(m, 1, y, t, tol);

	// y = alpha * A * x

	gemv_n(alpha, A, x, y);
	r = t * alpha;
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t + y * beta;
	gemv_n(A, x, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = alpha * A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t * alpha + y * beta;
	gemv_n(alpha, A, x, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);
}


TMN_CASE( gemv_n, refex )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;
	const T alpha = T(1.5);
	const T beta = T(1.2);
	T tol = tolvalue<T>::get();

	scoped_array<T> Abuf(LDim * n);

	ref_matrix_ex<T, M, N> A(Abuf.ptr_begin(), m, n, LDim);
	dense_matrix<T, N, 1> x(n, 1);
	dense_matrix<T, M, 1> y(m, 1, zeros<T>());
	dense_matrix<T, M, 1> r(m, 1, zeros<T>());
	dense_matrix<T, M, 1> t(m, 1, zeros<T>());

	fill_rand(A, T(-1), T(1));
	fill_rand(x, T(-1), T(1));

	// y = A * x

	gemv_n(A, x, y);
	naive_mtimes(A, x, t);
	ASSERT_MAT_APPROX(m, 1, y, t, tol);

	// y = alpha * A * x

	gemv_n(alpha, A, x, y);
	r = t * alpha;
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t + y * beta;
	gemv_n(A, x, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = alpha * A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t * alpha + y * beta;
	gemv_n(alpha, A, x, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);
}


TMN_CASE( gemv_n, mexpr )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;
	const T alpha = T(1.5);
	const T beta = T(1.2);
	T tol = tolvalue<T>::get();

	dense_matrix<T, M, N> A1(m, n);
	dense_matrix<T, M, N> A2(m, n);
	dense_matrix<T, N, 1> x1(n, 1);
	dense_matrix<T, N, 1> x2(n, 1);

	fill_rand(A1, T(-1), T(1));
	fill_rand(A2, T(-1), T(1));

	fill_rand(x1, T(-1), T(1));
	fill_rand(x2, T(-1), T(1));

	dense_matrix<T, M, N> A = A1 + A2;
	dense_matrix<T, N, 1> x = x1 + x2;
	dense_matrix<T, M, 1> y(m, 1, zeros<T>());
	dense_matrix<T, M, 1> r(m, 1, zeros<T>());
	dense_matrix<T, M, 1> t(m, 1, zeros<T>());

	// y = A * x

	gemv_n(A1 + A2, x1 + x2, y);
	naive_mtimes(A, x, t);
	ASSERT_MAT_APPROX(m, 1, y, t, tol);

	// y = alpha * A * x

	gemv_n(alpha, A1 + A2, x1 + x2, y);
	r = t * alpha;
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t + y * beta;
	gemv_n(A1 + A2, x1 + x2, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);

	// y = alpha * A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t * alpha + y * beta;
	gemv_n(alpha, A1 + A2, x1 + x2, beta, y);
	ASSERT_MAT_APPROX(m, 1, y, r, tol);
}


TMN_CASE( gemv_t, dense )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;
	const T alpha = T(1.5);
	const T beta = T(1.2);
	T tol = tolvalue<T>::get();

	dense_matrix<T, M, N> A(m, n);
	dense_matrix<T, M, 1> x(m, 1);
	dense_matrix<T, N, 1> y(n, 1, zeros<T>());
	dense_matrix<T, N, 1> r(n, 1, zeros<T>());
	dense_matrix<T, N, 1> t(n, 1, zeros<T>());

	fill_rand(A, T(-1), T(1));
	fill_rand(x, T(-1), T(1));

	dense_matrix<T, N, M> AT = A.trans();

	// y = A * x

	gemv_t(A, x, y);
	naive_mtimes(AT, x, t);
	ASSERT_MAT_APPROX(n, 1, y, t, tol);

	// y = alpha * A * x

	gemv_t(alpha, A, x, y);
	r = t * alpha;
	ASSERT_MAT_APPROX(n, 1, y, r, tol);

	// y = A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t + y * beta;
	gemv_t(A, x, beta, y);
	ASSERT_MAT_APPROX(n, 1, y, r, tol);

	// y = alpha * A * x + beta * y

	fill_rand(y, T(-1), T(1));
	r = t * alpha + y * beta;
	gemv_t(alpha, A, x, beta, y);
	ASSERT_MAT_APPROX(n, 1, y, r, tol);

}



BEGIN_TPACK( gemv_n_dense )
	ADD_TMN_CASE_3X3( gemv_n, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( gemv_n, dense, double, DM, DN )
END_TPACK

BEGIN_TPACK( gemv_n_refex )
	ADD_TMN_CASE_3X3( gemv_n, refex, float,  DM, DN )
	ADD_TMN_CASE_3X3( gemv_n, refex, double, DM, DN )
END_TPACK

BEGIN_TPACK( gemv_n_mexpr )
	ADD_TMN_CASE_3X3( gemv_n, mexpr, float,  DM, DN )
	ADD_TMN_CASE_3X3( gemv_n, mexpr, double, DM, DN )
END_TPACK

BEGIN_TPACK( gemv_t_dense )
	ADD_TMN_CASE_3X3( gemv_t, dense, float,  DM, DN )
	ADD_TMN_CASE_3X3( gemv_t, dense, double, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( gemv_n_dense )
	ADD_TPACK( gemv_n_refex )
	ADD_TPACK( gemv_n_mexpr )
	ADD_TPACK( gemv_t_dense )
END_MAIN_SUITE
