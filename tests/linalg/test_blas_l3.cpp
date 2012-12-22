/**
 * @file test_blas_l3.cpp
 *
 * @brief Unit testing of BLAS Level 3
 *
 * @author Dahua Lin
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l3.h>

using namespace lmat;
using namespace lmat::test;

const index_t DK = 5;

template<class SA, class SB, class SC, typename T, int M, int N>
void test_gemm_nn()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;
	index_t k = DK;

	typedef typename mat_host<SA, T, M, 0>::cmat_t amat_t;
	typedef typename mat_host<SB, T, 0, N>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, M, 0> a_host(m, k);
	mat_host<SB, T, 0, N> b_host(k, n);
	mat_host<SC, T, M, N> c_host(m, n);

	a_host.fill_rand();
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();


	safe_mm(T(1), a, 'n', b, 'n', T(0), c, r);
	blas::gemm(a, b, c);

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, a, 'n', b, 'n', beta, c, r);
	blas::gemm(alpha, a, b, beta, c);

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_gemm_nt()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;
	index_t k = DK;

	typedef typename mat_host<SA, T, M, 0>::cmat_t amat_t;
	typedef typename mat_host<SB, T, N, 0>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, M, 0> a_host(m, k);
	mat_host<SB, T, N, 0> b_host(n, k);
	mat_host<SC, T, M, N> c_host(m, n);

	a_host.fill_rand();
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), a, 'n', b, 't', T(0), c, r);
	blas::gemm(a, b, c, 'n', 't');

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, a, 'n', b, 't', beta, c, r);
	blas::gemm(alpha, a, b, beta, c, 'n', 't');

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_gemm_tn()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;
	index_t k = DK;

	typedef typename mat_host<SA, T, 0, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, 0, N>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, 0, M> a_host(k, m);
	mat_host<SB, T, 0, N> b_host(k, n);
	mat_host<SC, T, M, N> c_host(m, n);

	a_host.fill_rand();
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();


	safe_mm(T(1), a, 't', b, 'n', T(0), c, r);
	blas::gemm(a, b, c, 't', 'n');

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, a, 't', b, 'n', beta, c, r);
	blas::gemm(alpha, a, b, beta, c, 't', 'n');

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_gemm_tt()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;
	index_t k = DK;

	typedef typename mat_host<SA, T, 0, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, N, 0>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, 0, M> a_host(k, m);
	mat_host<SB, T, N, 0> b_host(n, k);
	mat_host<SC, T, M, N> c_host(m, n);

	a_host.fill_rand();
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), a, 't', b, 't', T(0), c, r);
	blas::gemm(a, b, c, 't', 't');

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, a, 't', b, 't', beta, c, r);
	blas::gemm(alpha, a, b, beta, c, 't', 't');

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}


#define DEF_BLAS3_CASES_P3( Name ) \
	TMN_CASE( mat_blas3, Name##_ccc ) { test_##Name<cont, cont, cont, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_ccb ) { test_##Name<cont, cont, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_cbc ) { test_##Name<cont, bloc, cont, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_cbb ) { test_##Name<cont, bloc, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_bcc ) { test_##Name<bloc, cont, cont, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_bcb ) { test_##Name<bloc, cont, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_bbc ) { test_##Name<bloc, bloc, cont, T, M, N>(); } \
	TMN_CASE( mat_blas3, Name##_bbb ) { test_##Name<bloc, bloc, bloc, T, M, N>(); }

#define ADD_BLAS3_CASES_P3( Name, T ) \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_ccc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_ccb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_cbc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_cbb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_bcc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_bcb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_bbc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas3, Name##_bbb, T, DM, DN );


// gemm

DEF_BLAS3_CASES_P3( gemm_nn )

BEGIN_TPACK( mat_gemm_nn_float )
	ADD_BLAS3_CASES_P3( gemm_nn, float )
END_TPACK

BEGIN_TPACK( mat_gemm_nn_double )
	ADD_BLAS3_CASES_P3( gemm_nn, double )
END_TPACK

DEF_BLAS3_CASES_P3( gemm_nt )

BEGIN_TPACK( mat_gemm_nt_float )
	ADD_BLAS3_CASES_P3( gemm_nt, float )
END_TPACK

BEGIN_TPACK( mat_gemm_nt_double )
	ADD_BLAS3_CASES_P3( gemm_nt, double )
END_TPACK

DEF_BLAS3_CASES_P3( gemm_tn )

BEGIN_TPACK( mat_gemm_tn_float )
	ADD_BLAS3_CASES_P3( gemm_tn, float )
END_TPACK

BEGIN_TPACK( mat_gemm_tn_double )
	ADD_BLAS3_CASES_P3( gemm_tn, double )
END_TPACK

DEF_BLAS3_CASES_P3( gemm_tt )

BEGIN_TPACK( mat_gemm_tt_float )
	ADD_BLAS3_CASES_P3( gemm_tt, float )
END_TPACK

BEGIN_TPACK( mat_gemm_tt_double )
	ADD_BLAS3_CASES_P3( gemm_tt, double )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_gemm_nn_float )
	ADD_TPACK( mat_gemm_nn_double )
	ADD_TPACK( mat_gemm_nt_float )
	ADD_TPACK( mat_gemm_nt_double )
	ADD_TPACK( mat_gemm_tn_float )
	ADD_TPACK( mat_gemm_tn_double )
	ADD_TPACK( mat_gemm_tt_float )
	ADD_TPACK( mat_gemm_tt_double )
END_MAIN_SUITE


