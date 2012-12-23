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
		dense_matrix<T> L(m, m);
		chol.get(L);
		blas::gemm(L, L, prod, 'n', 't');
	}
	else
	{
		dense_matrix<T> R(m, m);
		chol.get(R);
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



T_CASE( mat_chol, solve_l )
{
	test_chol_solve<T>('L');
}

T_CASE( mat_chol, solve_u )
{
	test_chol_solve<T>('U');
}

T_CASE( mat_chol, inv_l )
{
	test_chol_inv<T>('L');
}

T_CASE( mat_chol, inv_u )
{
	test_chol_inv<T>('U');
}

T_CASE( mat_chol, pdinv )
{
	test_pdinv<T>();
}

T_CASE( mat_chol, posv_l )
{
	test_posv<T>('L');
}

T_CASE( mat_chol, posv_u )
{
	test_posv<T>('U');
}


BEGIN_TPACK( mat_chol_solve )
	ADD_T_CASE( mat_chol, solve_l, float )
	ADD_T_CASE( mat_chol, solve_l, double )
	ADD_T_CASE( mat_chol, solve_u, float )
	ADD_T_CASE( mat_chol, solve_u, double )
END_TPACK


BEGIN_TPACK( mat_chol_inv )
	ADD_T_CASE( mat_chol, inv_l, float )
	ADD_T_CASE( mat_chol, inv_l, double )
	ADD_T_CASE( mat_chol, inv_u, float )
	ADD_T_CASE( mat_chol, inv_u, double )
END_TPACK

BEGIN_TPACK( mat_pdinv )
	ADD_T_CASE( mat_chol, pdinv, float )
	ADD_T_CASE( mat_chol, pdinv, double )
END_TPACK

BEGIN_TPACK( mat_posv )
	ADD_T_CASE( mat_chol, posv_l, float )
	ADD_T_CASE( mat_chol, posv_l, double )
	ADD_T_CASE( mat_chol, posv_u, float )
	ADD_T_CASE( mat_chol, posv_u, double )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_chol_solve )
	ADD_TPACK( mat_chol_inv )
	ADD_TPACK( mat_pdinv )
	ADD_TPACK( mat_posv )
END_MAIN_SUITE






