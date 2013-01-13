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

	dense_col<T, M> r(m);

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

	dense_row<T, N> r(n);

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


template<class SA, class SX, class SY, typename T, int N>
void test_symv_n()
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SX, T, N, 1>::cmat_t svec_t;
	typedef typename mat_host<SY, T, N, 1>::mat_t  dvec_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SX, T, N, 1> x_host(n, 1);
	mat_host<SY, T, N, 1> y_host(n, 1);

	mat_t amat = a_host.get_mat();
	fill_rand_sym(amat);
	x_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	svec_t x = x_host.get_cmat();
	dvec_t y = y_host.get_mat();

	dense_col<T, N> r(n);

	T tol = blas_default_tol<T>::get();

	safe_mv(T(1), a, 'n', x, T(0), y, r);
	blas::symv(a, x, y);

	ASSERT_VEC_APPROX(n, y, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mv(alpha, a, 'n', x, beta, y, r);
	blas::symv(alpha, a, x, beta, y);

	ASSERT_VEC_APPROX(n, y, r, tol);
}

template<class SA, class SX, class SY, typename T, int N>
void test_symv_t()
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SX, T, 1, N>::cmat_t svec_t;
	typedef typename mat_host<SY, T, 1, N>::mat_t  dvec_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SX, T, 1, N> x_host(1, n);
	mat_host<SY, T, 1, N> y_host(1, n);

	mat_t amat = a_host.get_mat();
	fill_rand_sym(amat);
	x_host.fill_rand();

	cmat_t a = a_host.get_cmat();
	svec_t x = x_host.get_cmat();
	dvec_t y = y_host.get_mat();

	dense_row<T, N> r(n);

	T tol = blas_default_tol<T>::get();

	safe_mv(T(1), a, 't', x, T(0), y, r);
	blas::symv(a, x, y, 'u');

	ASSERT_VEC_APPROX(n, y, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mv(alpha, a, 't', x, beta, y, r);
	blas::symv(alpha, a, x, beta, y, 'u');

	ASSERT_VEC_APPROX(n, y, r, tol);
}


template<class SA, class SY, typename T, int N>
void test_trmv_n(char uplo)
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SY, T, N, 1>::mat_t  dvec_t;

	T tol = blas_default_tol<T>::get();

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SY, T, N, 1> y_host(n, 1);

	mat_t amat = a_host.get_mat();
	cmat_t a = a_host.get_cmat();
	dvec_t y = y_host.get_mat();
	dense_col<T, N> r(n);

	fill_rand_tri(amat, uplo, false);
	y_host.fill_rand();

	safe_mv(T(1), a, 'n', y, T(0), y, r);
	blas::trmv(a, y, blas::trs(uplo, 'n'));

	ASSERT_VEC_APPROX(n, y, r, tol);

	fill_rand_tri(amat, uplo, true);
	y_host.fill_rand();

	safe_mv(T(1), a, 'n', y, T(0), y, r);
	blas::trmv(a, y, blas::trs(uplo, 'n', 'u'));

	ASSERT_VEC_APPROX(n, y, r, tol);
}

template<class SA, class SY, typename T, int N>
void test_trmv_t(char uplo)
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SY, T, 1, N>::mat_t  dvec_t;

	T tol = blas_default_tol<T>::get();

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SY, T, 1, N> y_host(1, n);

	mat_t amat = a_host.get_mat();
	cmat_t a = a_host.get_cmat();
	dvec_t y = y_host.get_mat();
	dense_row<double, N> r(n);

	fill_rand_tri(amat, uplo, false);
	y_host.fill_rand();

	safe_mv(T(1), a, 't', y, T(0), y, r);
	blas::trmv(a, y, blas::trs(uplo, 't'));

	ASSERT_VEC_APPROX(n, y, r, tol);

	fill_rand_tri(amat, uplo, true);
	y_host.fill_rand();

	safe_mv(T(1), a, 't', y, T(0), y, r);
	blas::trmv(a, y, blas::trs(uplo, 't', 'u'));

	ASSERT_VEC_APPROX(n, y, r, tol);
}


template<class SA, class SY, typename T, int N>
void test_trsv_n(char uplo)
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SY, T, N, 1>::mat_t  dvec_t;

	T tol = blas_default_tol<T>::get() * T(5);

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SY, T, N, 1> y_host(n, 1);

	mat_t amat = a_host.get_mat();
	cmat_t a = a_host.get_cmat();
	dvec_t y = y_host.get_mat();

	fill_rand_tri(amat, uplo, false);
	y_host.fill_rand();

	dense_col<T, N> p(n, zero());
	dense_col<T, N> r = y;

	blas::trsv(a, y, blas::trs(uplo, 'n'));
	safe_mv(T(1), a, 'n', y, T(0), y, p);

	ASSERT_VEC_APPROX(n, p, r, tol);

	fill_rand_tri(amat, uplo, true);
	y_host.fill_rand();

	dense_col<T, N> pu(n, zero());
	dense_col<T, N> ru = y;

	blas::trsv(a, y, blas::trs(uplo, 'n', 'u'));
	safe_mv(T(1), a, 'n', y, T(0), y, pu);

	ASSERT_VEC_APPROX(n, pu, ru, tol);
}


template<class SA, class SY, typename T, int N>
void test_trsv_t(char uplo)
{
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::mat_t mat_t;
	typedef typename mat_host<SA, T, N, N>::cmat_t cmat_t;
	typedef typename mat_host<SY, T, N, 1>::mat_t  dvec_t;

	T tol = blas_default_tol<T>::get() * T(5);

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SY, T, N, 1> y_host(n, 1);

	mat_t amat = a_host.get_mat();
	cmat_t a = a_host.get_cmat();
	dvec_t y = y_host.get_mat();

	fill_rand_tri(amat, uplo, false);
	y_host.fill_rand();

	dense_col<T, N> p(n, zero());
	dense_col<T, N> r = y;

	blas::trsv(a, y, blas::trs(uplo, 't'));
	safe_mv(T(1), a, 't', y, T(0), y, p);

	ASSERT_VEC_APPROX(n, p, r, tol);

	fill_rand_tri(amat, uplo, true);
	y_host.fill_rand();

	dense_col<T, N> pu(n, zero());
	dense_col<T, N> ru = y;

	blas::trsv(a, y, blas::trs(uplo, 't', 'u'));
	safe_mv(T(1), a, 't', y, T(0), y, pu);

	ASSERT_VEC_APPROX(n, pu, ru, tol);
}


template<class SA, class SX, class SY, typename T, int M, int N>
void test_ger()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, N>::mat_t mat_t;
	typedef typename mat_host<SX, T, M, 1>::cmat_t xvec_t;
	typedef typename mat_host<SY, T, 1, N>::cmat_t yvec_t;

	T tol = blas_default_tol<T>::get();

	mat_host<SA, T, M, N> a_host(m, n);
	mat_host<SX, T, M, 1> x_host(m, 1);
	mat_host<SY, T, 1, N> y_host(1, n);

	a_host.fill_rand();
	x_host.fill_rand();
	y_host.fill_rand();

	T alpha = T(3.2);

	mat_t a = a_host.get_mat();
	xvec_t x = x_host.get_cmat();
	yvec_t y = y_host.get_cmat();

	dense_matrix<T, M, N> r(a);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			r(i, j) += alpha * x[i] * y[j];
		}
	}

	blas::ger(a, x, y, alpha);

	ASSERT_MAT_APPROX( m, n, a, r, tol );

}



#define ADD_BLAS2_CASES_P3( Name, T ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_ccc, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_ccb, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_cbc, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_cbb, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bcc, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bcb, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bbc, T, DM, DN ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bbb, T, DM, DN )

#define ADD_BLAS2_CASES_P3S( Name, T ) \
	ADD_TN_CASE_3( mat_blas_##Name##_ccc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_ccb, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_cbc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_cbb, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bcc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bcb, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bbc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bbb, T, DN )

#define ADD_BLAS2_CASES_P2S( Name, T ) \
	ADD_TN_CASE_3( mat_blas_##Name##_cc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_cb, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bc, T, DN ) \
	ADD_TN_CASE_3( mat_blas_##Name##_bb, T, DN )


// gemv

TMN_CASE( mat_blas_gemv_n_ccc ) { test_gemv_n<cont, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_ccb ) { test_gemv_n<cont, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_cbc ) { test_gemv_n<cont, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_cbb ) { test_gemv_n<cont, bloc, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_bcc ) { test_gemv_n<bloc, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_bcb ) { test_gemv_n<bloc, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_bbc ) { test_gemv_n<bloc, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_n_bbb ) { test_gemv_n<bloc, bloc, bloc, T, M, N>(); }

TMN_CASE( mat_blas_gemv_t_ccc ) { test_gemv_t<cont, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_ccb ) { test_gemv_t<cont, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_cbc ) { test_gemv_t<cont, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_cbb ) { test_gemv_t<cont, bloc, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_bcc ) { test_gemv_t<bloc, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_bcb ) { test_gemv_t<bloc, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_bbc ) { test_gemv_t<bloc, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_gemv_t_bbb ) { test_gemv_t<bloc, bloc, bloc, T, M, N>(); }

AUTO_TPACK( mat_gemv_n )
{
	ADD_BLAS2_CASES_P3( gemv_n, float )
	ADD_BLAS2_CASES_P3( gemv_n, double )
}


AUTO_TPACK( mat_gemv_t )
{
	ADD_BLAS2_CASES_P3( gemv_t, float )
	ADD_BLAS2_CASES_P3( gemv_t, double )
}


// symv

TN_CASE( mat_blas_symv_n_ccc ) { test_symv_n<cont, cont, cont, T, N>(); }
TN_CASE( mat_blas_symv_n_ccb ) { test_symv_n<cont, cont, bloc, T, N>(); }
TN_CASE( mat_blas_symv_n_cbc ) { test_symv_n<cont, bloc, cont, T, N>(); }
TN_CASE( mat_blas_symv_n_cbb ) { test_symv_n<cont, bloc, bloc, T, N>(); }
TN_CASE( mat_blas_symv_n_bcc ) { test_symv_n<bloc, cont, cont, T, N>(); }
TN_CASE( mat_blas_symv_n_bcb ) { test_symv_n<bloc, cont, bloc, T, N>(); }
TN_CASE( mat_blas_symv_n_bbc ) { test_symv_n<bloc, bloc, cont, T, N>(); }
TN_CASE( mat_blas_symv_n_bbb ) { test_symv_n<bloc, bloc, bloc, T, N>(); }

TN_CASE( mat_blas_symv_t_ccc ) { test_symv_t<cont, cont, cont, T, N>(); }
TN_CASE( mat_blas_symv_t_ccb ) { test_symv_t<cont, cont, bloc, T, N>(); }
TN_CASE( mat_blas_symv_t_cbc ) { test_symv_t<cont, bloc, cont, T, N>(); }
TN_CASE( mat_blas_symv_t_cbb ) { test_symv_t<cont, bloc, bloc, T, N>(); }
TN_CASE( mat_blas_symv_t_bcc ) { test_symv_t<bloc, cont, cont, T, N>(); }
TN_CASE( mat_blas_symv_t_bcb ) { test_symv_t<bloc, cont, bloc, T, N>(); }
TN_CASE( mat_blas_symv_t_bbc ) { test_symv_t<bloc, bloc, cont, T, N>(); }
TN_CASE( mat_blas_symv_t_bbb ) { test_symv_t<bloc, bloc, bloc, T, N>(); }


AUTO_TPACK( mat_symv_n )
{
	ADD_BLAS2_CASES_P3S( symv_n, float )
	ADD_BLAS2_CASES_P3S( symv_n, double )
}

AUTO_TPACK( mat_symv_t )
{
	ADD_BLAS2_CASES_P3S( symv_t, float )
	ADD_BLAS2_CASES_P3S( symv_t, double )
}




// trmv

TN_CASE( mat_blas_trmv_n_cc ) { test_trmv_n<cont, cont, T, N>('l'); test_trmv_n<cont, cont, T, N>('u'); }
TN_CASE( mat_blas_trmv_n_cb ) { test_trmv_n<cont, bloc, T, N>('l'); test_trmv_n<cont, bloc, T, N>('u'); }
TN_CASE( mat_blas_trmv_n_bc ) { test_trmv_n<bloc, cont, T, N>('l'); test_trmv_n<bloc, cont, T, N>('u'); }
TN_CASE( mat_blas_trmv_n_bb ) { test_trmv_n<bloc, bloc, T, N>('l'); test_trmv_n<bloc, bloc, T, N>('u'); }

TN_CASE( mat_blas_trmv_t_cc ) { test_trmv_t<cont, cont, T, N>('l'); test_trmv_t<cont, cont, T, N>('u'); }
TN_CASE( mat_blas_trmv_t_cb ) { test_trmv_t<cont, bloc, T, N>('l'); test_trmv_t<cont, bloc, T, N>('u'); }
TN_CASE( mat_blas_trmv_t_bc ) { test_trmv_t<bloc, cont, T, N>('l'); test_trmv_t<bloc, cont, T, N>('u'); }
TN_CASE( mat_blas_trmv_t_bb ) { test_trmv_t<bloc, bloc, T, N>('l'); test_trmv_t<bloc, bloc, T, N>('u'); }

AUTO_TPACK( mat_trmv_n )
{
	ADD_BLAS2_CASES_P2S( trmv_n, float )
	ADD_BLAS2_CASES_P2S( trmv_n, double )
}


AUTO_TPACK( mat_trmv_t )
{
	ADD_BLAS2_CASES_P2S( trmv_t, float )
	ADD_BLAS2_CASES_P2S( trmv_t, double )
}


// trsv

TN_CASE( mat_blas_trsv_n_cc ) { test_trsv_n<cont, cont, T, N>('l'); test_trsv_n<cont, cont, T, N>('u'); }
TN_CASE( mat_blas_trsv_n_cb ) { test_trsv_n<cont, bloc, T, N>('l'); test_trsv_n<cont, bloc, T, N>('u'); }
TN_CASE( mat_blas_trsv_n_bc ) { test_trsv_n<bloc, cont, T, N>('l'); test_trsv_n<bloc, cont, T, N>('u'); }
TN_CASE( mat_blas_trsv_n_bb ) { test_trsv_n<bloc, bloc, T, N>('l'); test_trsv_n<bloc, bloc, T, N>('u'); }

TN_CASE( mat_blas_trsv_t_cc ) { test_trsv_t<cont, cont, T, N>('l'); test_trsv_t<cont, cont, T, N>('u'); }
TN_CASE( mat_blas_trsv_t_cb ) { test_trsv_t<cont, bloc, T, N>('l'); test_trsv_t<cont, bloc, T, N>('u'); }
TN_CASE( mat_blas_trsv_t_bc ) { test_trsv_t<bloc, cont, T, N>('l'); test_trsv_t<bloc, cont, T, N>('u'); }
TN_CASE( mat_blas_trsv_t_bb ) { test_trsv_t<bloc, bloc, T, N>('l'); test_trsv_t<bloc, bloc, T, N>('u'); }

AUTO_TPACK( mat_trsv_n )
{
	ADD_BLAS2_CASES_P2S( trsv_n, float )
	ADD_BLAS2_CASES_P2S( trsv_n, double )
}

AUTO_TPACK( mat_trsv_t )
{
	ADD_BLAS2_CASES_P2S( trsv_t, float )
	ADD_BLAS2_CASES_P2S( trsv_t, double )
}



// ger

TMN_CASE( mat_blas_ger_ccc ) { test_ger<cont, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_ger_ccb ) { test_ger<cont, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_ger_cbc ) { test_ger<cont, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_ger_cbb ) { test_ger<cont, bloc, bloc, T, M, N>(); }
TMN_CASE( mat_blas_ger_bcc ) { test_ger<bloc, cont, cont, T, M, N>(); }
TMN_CASE( mat_blas_ger_bcb ) { test_ger<bloc, cont, bloc, T, M, N>(); }
TMN_CASE( mat_blas_ger_bbc ) { test_ger<bloc, bloc, cont, T, M, N>(); }
TMN_CASE( mat_blas_ger_bbb ) { test_ger<bloc, bloc, bloc, T, M, N>(); }

AUTO_TPACK( mat_ger )
{
	ADD_BLAS2_CASES_P3( ger, float )
	ADD_BLAS2_CASES_P3( ger, double )
}



