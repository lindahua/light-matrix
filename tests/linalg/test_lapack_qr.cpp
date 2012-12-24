/**
 * @file test_lapack_qr.cpp
 *
 * @brief Unit testing of QR factorization
 *
 * @author Dahua Lin
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l3.h>
#include <light_mat/linalg/lapack_qr.h>

using namespace lmat;
using namespace lmat::test;
using lmat::lapack::qr_fac;


template<typename T>
void test_qr_fac(index_t m, index_t n)
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::mat_t mat_t;
	typedef typename host_t::cmat_t cmat_t;

	host_t a_host(m, n);
	mat_t a = a_host.get_mat();

	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);

	randunif(T(0), T(1));

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) a(i, j) = randunif(T(-2.0), T(2.0));
	}

	qr_fac<T> qr(a);

	dense_matrix<T> r0(m, n, zero());
	qr.getr(r0);

	dense_matrix<T> q0(m, m, zero());
	qr.getq(q0);

	dense_matrix<T> e0(m, m); fill_eye(e0);
	dense_matrix<T> g0(m, m); blas::gemm(q0, q0, g0, 'T', 'N');
	ASSERT_MAT_APPROX(m, m, g0, e0, tol);

	dense_matrix<T> p0(m, n, zero());
	blas::gemm(q0, r0, p0);
	ASSERT_MAT_APPROX(m, n, p0, a, tol);

	if (m > n)
	{
		dense_matrix<T> r0s(n, n, zero());
		qr.getr(r0s);
		dense_matrix<T> q0s(m, n, zero());
		qr.getq(q0s);

		ASSERT_MAT_EQ( m, n, q0, q0s );

		dense_matrix<T> p0s(m, n, zero());
		blas::gemm(q0s, r0s, p0s);
		ASSERT_MAT_APPROX(m, n, p0s, a, tol);
	}

	index_t n2 = math::min(m, n) - 2;
	dense_matrix<T> q2(m, n2);
	qr.getq(q2);
	ASSERT_MAT_EQ( m, n2, q0, q2 );
}


T_CASE( mat_qr, fac_eq )
{
	test_qr_fac<T>(5, 5);
}

T_CASE( mat_qr, fac_gt )
{
	test_qr_fac<T>(8, 5);
}

T_CASE( mat_qr, fac_lt )
{
	test_qr_fac<T>(5, 8);
}


BEGIN_TPACK( mat_qr_fac )
	ADD_T_CASE( mat_qr, fac_eq, float )
	ADD_T_CASE( mat_qr, fac_gt, float )
	ADD_T_CASE( mat_qr, fac_lt, float )
	ADD_T_CASE( mat_qr, fac_eq, double )
	ADD_T_CASE( mat_qr, fac_gt, double )
	ADD_T_CASE( mat_qr, fac_lt, double )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_qr_fac )
END_MAIN_SUITE
