/**
 * @file test_blas_l2.cpp
 *
 * @brief Unit testing of BLAS Level 2
 *
 * @author Dahua Lin
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l2.h>

using namespace lmat;
using namespace lmat::test;


template<class SA, class SX, class SY, typename T, int M, int N>
void test_gemv_n()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, N>::cmat_t cmat_t;
	typedef typename mat_host<SX, T, N, 1>::cmat_t svec_t;
	typedef typename mat_host<SY, T, M, 1>::mat_t  dvec_t;

	mat_host<SA, T, M, N> a_host(m, n);
	mat_host<SX, T, N, 1> x_host(n, 1);
	mat_host<SY, T, M, 1> y_host(m, 1);

	a_host.fill_rand();
	x_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	svec_t x = x_host.get_cmat();
	dvec_t y = y_host.get_mat();

	dense_col<double, M> r(m);

	T tol = blas_default_tol<T>::get();

	safe_mv(T(1), a, 'n', x, T(0), y, r);
	blas::gemv(a, x, y);

	ASSERT_VEC_APPROX(m, y, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mv(alpha, a, 'n', x, beta, y, r);
	blas::gemv(alpha, a, x, beta, y);

	ASSERT_VEC_APPROX(m, y, r, tol);
}

template<class SA, class SX, class SY, typename T, int M, int N>
void test_gemv_t()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, N>::cmat_t cmat_t;
	typedef typename mat_host<SX, T, 1, M>::cmat_t svec_t;
	typedef typename mat_host<SY, T, 1, N>::mat_t  dvec_t;

	mat_host<SA, T, M, N> a_host(m, n);
	mat_host<SX, T, 1, M> x_host(1, m);
	mat_host<SY, T, 1, N> y_host(1, n);

	a_host.fill_rand();
	x_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	svec_t x = x_host.get_cmat();
	dvec_t y = y_host.get_mat();

	dense_row<double, N> r(n);

	T tol = blas_default_tol<T>::get();

	safe_mv(T(1), a, 't', x, T(0), y, r);
	blas::gemv(a, x, y, 't');

	ASSERT_VEC_APPROX(n, y, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mv(alpha, a, 't', x, beta, y, r);
	blas::gemv(alpha, a, x, beta, y, 't');

	ASSERT_VEC_APPROX(n, y, r, tol);
}



TMN_CASE( mat_blas2, gemv_n_ccc ) { test_gemv_n<cont, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_ccb ) { test_gemv_n<cont, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_cbc ) { test_gemv_n<cont, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_cbb ) { test_gemv_n<cont, bloc, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_bcc ) { test_gemv_n<bloc, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_bcb ) { test_gemv_n<bloc, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_bbc ) { test_gemv_n<bloc, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_n_bbb ) { test_gemv_n<bloc, bloc, bloc, T, M, N>(); }

TMN_CASE( mat_blas2, gemv_t_ccc ) { test_gemv_t<cont, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_ccb ) { test_gemv_t<cont, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_cbc ) { test_gemv_t<cont, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_cbb ) { test_gemv_t<cont, bloc, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_bcc ) { test_gemv_t<bloc, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_bcb ) { test_gemv_t<bloc, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_bbc ) { test_gemv_t<bloc, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas2, gemv_t_bbb ) { test_gemv_t<bloc, bloc, bloc, T, M, N>(); }

BEGIN_TPACK( mat_gemv_n )
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_ccc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_ccb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_cbc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_cbb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_bcc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_bcb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_bbc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_n_bbb, double, DM, DN );
END_TPACK

BEGIN_TPACK( mat_gemv_t )
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_ccc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_ccb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_cbc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_cbb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_bcc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_bcb, double, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_bbc, float, DM, DN );
	ADD_TMN_CASE_3X3( mat_blas2, gemv_t_bbb, double, DM, DN );
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_gemv_n )
	ADD_TPACK( mat_gemv_t )
END_MAIN_SUITE

