/**
 * @file test_lapack_chol.cpp
 *
 * @brief Unit testing of Cholesky decomposition
 *
 * @author Dahua Lin
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l3.h>
#include <light_mat/linalg/lapack_chol.h>


using namespace lmat;
using namespace lmat::test;

using lmat::lapack::chol_fac;

template<typename T>
void test_chol_solve( char uplo )
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	index_t m = DM;
	index_t n = DN;

	host_t a_host(m, m);
	host_t b_host(m, n);
	host_t x_host(m, n);

	mat_t amat = a_host.get_mat();
	fill_rand_pdm(amat);
	b_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	cmat_t b = b_host.get_cmat();
	mat_t x = x_host.get_mat();

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);

	chol_fac<T> chol(a, uplo);
	ASSERT_EQ( chol.dim(), m );
	ASSERT_EQ( chol.is_lower(), uplo == 'L' );
	ASSERT_EQ( chol.is_upper(), uplo == 'U' );

	dense_matrix<T> prod(m, m);

	if (chol.is_lower())
	{
		dense_matrix<T> L;
		chol.get(L);
		ASSERT_EQ(L.nrows(), m);
		ASSERT_EQ(L.ncolumns(), m);

		blas::gemm(L, L, prod, 'n', 't');
	}
	else
	{
		dense_matrix<T> R(m, m);
		chol.get(R);
		ASSERT_EQ(R.nrows(), m);
		ASSERT_EQ(R.ncolumns(), m);

		blas::gemm(R, R, prod, 't', 'n');
	}

	ASSERT_MAT_APPROX( m, m, prod, a, tol );

	chol.solve(b, x);

	dense_matrix<T> r(m, n);
	blas::gemm(a, x, r);

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}


template<typename T>
void test_chol_inv( char uplo )
{
	typedef mat_host<cont, T, 0, 0> host_t;
	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	index_t m = DM;

	host_t a_host(m, m);
	mat_t amat = a_host.get_mat();
	fill_rand_pdm(amat);

	cmat_t a = a_host.get_cmat();

	dense_matrix<T> b(m, m);
	chol_fac<T>::inv(a, b, uplo);

	dense_matrix<T> r(m, m);
	dense_matrix<T> e(m, m);

	blas::gemm(a, b, r);
	fill_eye(e);

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);
	ASSERT_MAT_APPROX(m, m, r, e, tol);
}


template<typename T>
void test_pdinv()
{
	typedef mat_host<cont, T, 0, 0> host_t;
	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	index_t m = DM;

	host_t a_host(m, m);
	mat_t amat = a_host.get_mat();
	fill_rand_pdm(amat);

	cmat_t a = a_host.get_cmat();

	dense_matrix<T> b = pdinv(a);

	ASSERT_EQ( b.nrows(), m );
	ASSERT_EQ( b.ncolumns(), m );

	dense_matrix<T> r(m, m);
	dense_matrix<T> e(m, m);

	blas::gemm(a, b, r);
	fill_eye(e);

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);
	ASSERT_MAT_APPROX(m, m, r, e, tol);
}


template<typename T>
void test_posv( char uplo )
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	index_t m = DM;
	index_t n = DN;

	host_t a_host(m, m);
	host_t b_host(m, n);

	mat_t amat = a_host.get_mat();
	fill_rand_pdm(amat);
	b_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	cmat_t b = b_host.get_cmat();

	dense_matrix<T> a_(a);
	dense_matrix<T> x(b);
	lapack::posv(a_, x, uplo);

	dense_matrix<T> r(m, n);
	blas::gemm(a, x, r);

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}

template<typename T>
void test_pddet(const dense_matrix<T>& a, const T& vdet)
{
	const T lddet = math::log(vdet);

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10) * vdet;

	ASSERT_APPROX( pddet(a), vdet, tol );
	ASSERT_APPROX( pdlogdet(a), lddet, tol );
}


T_CASE( mat_chol_solve_l )
{
	test_chol_solve<T>('L');
}

T_CASE( mat_chol_solve_u )
{
	test_chol_solve<T>('U');
}

T_CASE( mat_chol_inv_l )
{
	test_chol_inv<T>('L');
}

T_CASE( mat_chol_inv_u )
{
	test_chol_inv<T>('U');
}

T_CASE( mat_chol_pdinv )
{
	test_pdinv<T>();
}

T_CASE( mat_posv_l )
{
	test_posv<T>('L');
}

T_CASE( mat_posv_u )
{
	test_posv<T>('U');
}


T_CASE( mat_pddet_1 )
{
	T v = T(12.3);
	dense_matrix<T> a(1, 1);
	a[0] = v;

	test_pddet(a, v);
}

T_CASE( mat_pddet_2 )
{
	int isrc[4] = { 3, -2, -2, 4 };
	T src[4];
	for (int i = 0; i < 4; ++i) src[i] = T(isrc[i]);

	dense_matrix<T> a(2, 2, copy_from(src));
	T v = T(8);

	test_pddet(a, v);
}

T_CASE( mat_pddet_3 )
{
	int isrc[9] = { 8, -2, -1, -2, 5, -3, -1, -3, 6 };
	T src[9];
	for (int i = 0; i < 9; ++i) src[i] = T(isrc[i]);

	dense_matrix<T> a(3, 3, copy_from(src));
	T v = T(127);

	test_pddet(a, v);
}


T_CASE( mat_pddet_5 )
{
	int isrc[25] = {
		 8,  0, -2,  1, -6,
		 0,  6,  1,  1,  2,
		-2,  1,  6,  1,  5,
		 1,  1,  1,  8,  1,
		-6,  2,  5,  1, 18 };

	T src[25];
	for (int i = 0; i < 25; ++i) src[i] = T(isrc[i]);

	dense_matrix<T> a(5, 5, copy_from(src));
	T v = T(20844);

	test_pddet(a, v);
}



AUTO_TPACK( mat_chol_solve )
{
	ADD_T_CASE( mat_chol_solve_l, float )
	ADD_T_CASE( mat_chol_solve_l, double )
	ADD_T_CASE( mat_chol_solve_u, float )
	ADD_T_CASE( mat_chol_solve_u, double )
}


AUTO_TPACK( mat_chol_inv )
{
	ADD_T_CASE( mat_chol_inv_l, float )
	ADD_T_CASE( mat_chol_inv_l, double )
	ADD_T_CASE( mat_chol_inv_u, float )
	ADD_T_CASE( mat_chol_inv_u, double )
}

AUTO_TPACK( mat_pdinv )
{
	ADD_T_CASE( mat_chol_pdinv, float )
	ADD_T_CASE( mat_chol_pdinv, double )
}

AUTO_TPACK( mat_posv )
{
	ADD_T_CASE( mat_posv_l, float )
	ADD_T_CASE( mat_posv_l, double )
	ADD_T_CASE( mat_posv_u, float )
	ADD_T_CASE( mat_posv_u, double )
}

AUTO_TPACK( mat_pddev )
{
	ADD_T_CASE( mat_pddet_1, float )
	ADD_T_CASE( mat_pddet_1, double )
	ADD_T_CASE( mat_pddet_2, float )
	ADD_T_CASE( mat_pddet_2, double )
	ADD_T_CASE( mat_pddet_3, float )
	ADD_T_CASE( mat_pddet_3, double )
	ADD_T_CASE( mat_pddet_5, float )
	ADD_T_CASE( mat_pddet_5, double )
}


