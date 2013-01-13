/**
 * @file test_blas_l1.cpp
 *
 * Unit testing of BLAS Level 1
 * 
 * @author Dahua Lin 
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l1.h>

using namespace lmat;
using namespace lmat::test;


template<typename S, typename T, int M, int N>
void test_asum()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::cmat_t smat_t;

	mat_host<S, T, M, N> x_host(m, n);
	x_host.fill_rand();

	smat_t x = x_host.get_cmat();
	index_t incx = safe_get_vec_intv(x);

	T r0 = safe_asum(m * n, x.ptr_data(), incx);
	T r = blas::asum(x);

	ASSERT_APPROX(r, r0, tol);
}


template<typename S, typename T, int M, int N>
void test_axpy()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::cmat_t smat_t;
	typedef typename mat_host<S, T, M, N>::mat_t dmat_t;

	mat_host<S, T, M, N> x_host(m, n);
	mat_host<S, T, M, N> y_host(m, n);
	x_host.fill_rand();
	y_host.fill_rand();

	T a = T(2.5);

	smat_t x = x_host.get_cmat();
	dmat_t y = y_host.get_mat();
	index_t incx = safe_get_vec_intv(x);
	index_t incy = safe_get_vec_intv(y);

	dense_matrix<T, M, N> r(m, n);
	safe_axpy(m * n, a, x.ptr_data(), incx, y.ptr_data(), incy, r.ptr_data());

	blas::axpy(a, x, y);

	ASSERT_MAT_APPROX(m, n, y, r, tol);
}


template<typename S, typename T, int M, int N>
void test_nrm2()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::cmat_t smat_t;

	mat_host<S, T, M, N> x_host(m, n);
	x_host.fill_rand();

	smat_t x = x_host.get_cmat();
	index_t incx = safe_get_vec_intv(x);

	T r0 = safe_nrm2(m * n, x.ptr_data(), incx);
	T r = blas::nrm2(x);

	ASSERT_APPROX(r, r0, tol);
}


template<typename S, typename T, int M, int N>
void test_dot()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::cmat_t smat_t;

	mat_host<S, T, M, N> x_host(m, n);
	mat_host<S, T, M, N> y_host(m, n);
	x_host.fill_rand();
	y_host.fill_rand();

	smat_t x = x_host.get_cmat();
	smat_t y = y_host.get_cmat();
	index_t incx = safe_get_vec_intv(x);
	index_t incy = safe_get_vec_intv(y);

	T r0 = safe_dot(m * n, x.ptr_data(), incx, y.ptr_data(), incy);
	T r = blas::dot(x, y);

	ASSERT_APPROX(r, r0, tol);
}


template<typename S, typename T, int M, int N>
void test_rot()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::mat_t dmat_t;

	mat_host<S, T, M, N> x_host(m, n);
	mat_host<S, T, M, N> y_host(m, n);
	x_host.fill_rand();
	y_host.fill_rand();

	T c = T(2.4);
	T s = T(1.6);

	dmat_t x = x_host.get_mat();
	dmat_t y = y_host.get_mat();
	index_t incx = safe_get_vec_intv(x);
	index_t incy = safe_get_vec_intv(y);

	dense_matrix<T, M, N> rx(m, n);
	dense_matrix<T, M, N> ry(m, n);
	safe_rot(m * n, c, s, x.ptr_data(), incx, y.ptr_data(), incy, rx.ptr_data(), ry.ptr_data());

	blas::rot(x, y, c, s);

	ASSERT_MAT_APPROX(m, n, x, rx, tol);
	ASSERT_MAT_APPROX(m, n, y, ry, tol);
}


template<typename S, typename T, int M, int N>
void test_scal()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	T tol = blas_default_tol<T>::get();

	typedef typename mat_host<S, T, M, N>::mat_t dmat_t;

	mat_host<S, T, M, N> x_host(m, n);
	mat_host<S, T, M, N> y_host(m, n);
	x_host.fill_rand();
	y_host.fill_rand();

	T a = T(3.2);

	dmat_t x = x_host.get_mat();
	index_t incx = safe_get_vec_intv(x);

	dense_matrix<T, M, N> r(m, n);
	safe_scal(m * n, a, x.ptr_data(), incx, r.ptr_data());

	blas::scal(x, a);

	ASSERT_MAT_APPROX(m, n, x, r, tol);
}




#define DEFINE_BLAS_L1_CASE( Name ) \
	TMN_CASE( mat_blas_##Name##_cont ) { test_##Name<cont, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_bloc ) { test_##Name<bloc, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_grid ) { test_##Name<grid, T, M, N>(); } \
	AUTO_TPACK( mat_blas_##Name ) { \
		ADD_TMN_CASE_3X3( mat_blas_##Name##_cont, float, DM, DN ) \
		ADD_TMN_CASE_3X3( mat_blas_##Name##_cont, double, DM, DN ) \
		ADD_TMN_CASE( mat_blas_##Name##_bloc, float, 1, 0 ) \
		ADD_TMN_CASE( mat_blas_##Name##_bloc, double, 1, DN ) \
		ADD_TMN_CASE( mat_blas_##Name##_grid, float, 0, 1 ) \
		ADD_TMN_CASE( mat_blas_##Name##_grid, double, DM, 1 ) \
	}

DEFINE_BLAS_L1_CASE( asum )
DEFINE_BLAS_L1_CASE( axpy )
DEFINE_BLAS_L1_CASE( nrm2 )
DEFINE_BLAS_L1_CASE( dot )
DEFINE_BLAS_L1_CASE( rot )
DEFINE_BLAS_L1_CASE( scal )

