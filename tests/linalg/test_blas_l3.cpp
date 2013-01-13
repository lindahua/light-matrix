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


template<class SA, class SB, class SC, typename T, int M, int N>
void test_symm_l()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, M, M> a_host(m, m);
	mat_host<SB, T, M, N> b_host(m, n);
	mat_host<SC, T, M, N> c_host(m, n);

	typename mat_host<SA, T, M, M>::mat_t amat = a_host.get_mat();
	fill_rand_sym(amat);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), a, 'n', b, 'n', T(0), c, r);
	blas::symm(a, b, c);

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, a, 'n', b, 'n', beta, c, r);
	blas::symm(alpha, a, b, beta, c);

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_symm_r()
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::cmat_t bmat_t;
	typedef typename mat_host<SC, T, M, N>::mat_t  cmat_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SB, T, M, N> b_host(m, n);
	mat_host<SC, T, M, N> c_host(m, n);

	typename mat_host<SA, T, N, N>::mat_t amat = a_host.get_mat();
	fill_rand_sym(amat);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_cmat();
	cmat_t c = c_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), b, 'n', a, 'n', T(0), c, r);
	blas::symm(a, b, c, 'r');

	ASSERT_MAT_APPROX(m, n, c, r, tol);

	T alpha = T(2.5);
	T beta = T(1.6);

	safe_mm(alpha, b, 'n', a, 'n', beta, c, r);
	blas::symm(alpha, a, b, beta, c, 'r');

	ASSERT_MAT_APPROX(m, n, c, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_trmm_ln(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, M, M> a_host(m, m);
	mat_host<SB, T, M, N> b_host(m, n);

	typename mat_host<SA, T, M, M>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), a, 'n', b, 'n', T(0), b, r);
	blas::trmm(a, b, blas::trs(uplo));

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_trmm_lt(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, M, M> a_host(m, m);
	mat_host<SB, T, M, N> b_host(m, n);

	typename mat_host<SA, T, M, M>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), a, 't', b, 'n', T(0), b, r);
	blas::trmm(a, b, blas::trs(uplo, 't'));

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_trmm_rn(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SB, T, M, N> b_host(m, n);

	typename mat_host<SA, T, N, N>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), b, 'n', a, 'n', T(0), b, r);
	blas::trmm(a, b, blas::trs(uplo), 'r');

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_trmm_rt(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SB, T, M, N> b_host(m, n);

	typename mat_host<SA, T, N, N>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r(m, n);
	T tol = blas_default_tol<T>::get();

	safe_mm(T(1), b, 'n', a, 't', T(0), b, r);
	blas::trmm(a, b, blas::trs(uplo, 't'), 'r');

	ASSERT_MAT_APPROX(m, n, b, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_trsm_ln(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, M, M> a_host(m, m);
	mat_host<SB, T, M, N> b_host(m, n);

	T tol = blas_default_tol<T>::get() * T(5);

	typename mat_host<SA, T, M, M>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r = b;
	dense_matrix<T> p(m, n, zero());

	blas::trsm(a, b, blas::trs(uplo));
	safe_mm(T(1), a, 'n', b, 'n', T(0), b, p);

	ASSERT_MAT_APPROX(m, n, p, r, tol);
}


template<class SA, class SB, class SC, typename T, int M, int N>
void test_trsm_lt(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, M, M>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, M, M> a_host(m, m);
	mat_host<SB, T, M, N> b_host(m, n);

	T tol = blas_default_tol<T>::get() * T(5);

	typename mat_host<SA, T, M, M>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r = b;
	dense_matrix<T> p(m, n, zero());

	blas::trsm(a, b, blas::trs(uplo, 't'));
	safe_mm(T(1), a, 't', b, 'n', T(0), b, p);

	ASSERT_MAT_APPROX(m, n, p, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_trsm_rn(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SB, T, M, N> b_host(m, n);

	T tol = blas_default_tol<T>::get() * T(5);

	typename mat_host<SA, T, N, N>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r = b;
	dense_matrix<T> p(m, n, zero());

	blas::trsm(a, b, blas::trs(uplo), 'r');
	safe_mm(T(1), b, 'n', a, 'n', T(0), b, p);

	ASSERT_MAT_APPROX(m, n, p, r, tol);
}

template<class SA, class SB, class SC, typename T, int M, int N>
void test_trsm_rt(char uplo)
{
	index_t m = M == 0 ? DM : M;
	index_t n = N == 0 ? DN : N;

	typedef typename mat_host<SA, T, N, N>::cmat_t amat_t;
	typedef typename mat_host<SB, T, M, N>::mat_t bmat_t;

	mat_host<SA, T, N, N> a_host(n, n);
	mat_host<SB, T, M, N> b_host(m, n);

	T tol = blas_default_tol<T>::get() * T(5);

	typename mat_host<SA, T, N, N>::mat_t amat = a_host.get_mat();
	fill_rand_tri(amat, uplo);
	b_host.fill_rand();

	amat_t a = a_host.get_cmat();
	bmat_t b = b_host.get_mat();

	dense_matrix<T> r = b;
	dense_matrix<T> p(m, n, zero());

	blas::trsm(a, b, blas::trs(uplo, 't'), 'r');
	safe_mm(T(1), b, 'n', a, 't', T(0), b, p);

	ASSERT_MAT_APPROX(m, n, p, r, tol);
}




#define DEF_BLAS3_CASES_P3( Name ) \
	TMN_CASE( mat_blas_##Name##_ccc ) { test_##Name<cont, cont, cont, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_ccb ) { test_##Name<cont, cont, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_cbc ) { test_##Name<cont, bloc, cont, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_cbb ) { test_##Name<cont, bloc, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_bcc ) { test_##Name<bloc, cont, cont, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_bcb ) { test_##Name<bloc, cont, bloc, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_bbc ) { test_##Name<bloc, bloc, cont, T, M, N>(); } \
	TMN_CASE( mat_blas_##Name##_bbb ) { test_##Name<bloc, bloc, bloc, T, M, N>(); }

#define DEF_BLAS3_CASES_P3_TR( Name ) \
	TMN_CASE( mat_blas_##Name##_ccc ) { test_##Name<cont, cont, cont, T, M, N>('l'); test_##Name<cont, cont, cont, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_ccb ) { test_##Name<cont, cont, bloc, T, M, N>('l'); test_##Name<cont, cont, bloc, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_cbc ) { test_##Name<cont, bloc, cont, T, M, N>('l'); test_##Name<cont, bloc, cont, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_cbb ) { test_##Name<cont, bloc, bloc, T, M, N>('l'); test_##Name<cont, bloc, bloc, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_bcc ) { test_##Name<bloc, cont, cont, T, M, N>('l'); test_##Name<bloc, cont, cont, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_bcb ) { test_##Name<bloc, cont, bloc, T, M, N>('l'); test_##Name<bloc, cont, bloc, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_bbc ) { test_##Name<bloc, bloc, cont, T, M, N>('l'); test_##Name<bloc, bloc, cont, T, M, N>('u'); } \
	TMN_CASE( mat_blas_##Name##_bbb ) { test_##Name<bloc, bloc, bloc, T, M, N>('l'); test_##Name<bloc, bloc, bloc, T, M, N>('u'); }

#define ADD_BLAS3_CASES_P3( Name, T ) \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_ccc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_ccb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_cbc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_cbb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bcc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bcb, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bbc, T, DM, DN ); \
	ADD_TMN_CASE_3X3( mat_blas_##Name##_bbb, T, DM, DN );


// gemm

DEF_BLAS3_CASES_P3( gemm_nn )

AUTO_TPACK( mat_gemm_nn )
{
	ADD_BLAS3_CASES_P3( gemm_nn, float )
	ADD_BLAS3_CASES_P3( gemm_nn, double )
}

DEF_BLAS3_CASES_P3( gemm_nt )

AUTO_TPACK( mat_gemm_nt )
{
	ADD_BLAS3_CASES_P3( gemm_nt, float )
	ADD_BLAS3_CASES_P3( gemm_nt, double )
}

DEF_BLAS3_CASES_P3( gemm_tn )

AUTO_TPACK( mat_gemm_tn )
{
	ADD_BLAS3_CASES_P3( gemm_tn, float )
	ADD_BLAS3_CASES_P3( gemm_tn, double )
}

DEF_BLAS3_CASES_P3( gemm_tt )

AUTO_TPACK( mat_gemm_tt )
{
	ADD_BLAS3_CASES_P3( gemm_tt, float )
	ADD_BLAS3_CASES_P3( gemm_tt, double )
}



// symm

DEF_BLAS3_CASES_P3( symm_l )

AUTO_TPACK( mat_symm_l )
{
	ADD_BLAS3_CASES_P3( symm_l, float )
	ADD_BLAS3_CASES_P3( symm_l, double )
}

DEF_BLAS3_CASES_P3( symm_r )

AUTO_TPACK( mat_symm_r )
{
	ADD_BLAS3_CASES_P3( symm_r, float )
	ADD_BLAS3_CASES_P3( symm_r, double )
}




// trmm

DEF_BLAS3_CASES_P3_TR( trmm_ln )

AUTO_TPACK( mat_trmm_ln )
{
	ADD_BLAS3_CASES_P3( trmm_ln, float )
	ADD_BLAS3_CASES_P3( trmm_ln, double )
}


DEF_BLAS3_CASES_P3_TR( trmm_lt )

AUTO_TPACK( mat_trmm_lt )
{
	ADD_BLAS3_CASES_P3( trmm_lt, float )
	ADD_BLAS3_CASES_P3( trmm_lt, double )
}


DEF_BLAS3_CASES_P3_TR( trmm_rn )

AUTO_TPACK( mat_trmm_rn )
{
	ADD_BLAS3_CASES_P3( trmm_rn, float )
	ADD_BLAS3_CASES_P3( trmm_rn, double )
}

DEF_BLAS3_CASES_P3_TR( trmm_rt )

AUTO_TPACK( mat_trmm_rt )
{
	ADD_BLAS3_CASES_P3( trmm_rt, float )
	ADD_BLAS3_CASES_P3( trmm_rt, double )
}


// tsmm

DEF_BLAS3_CASES_P3_TR( trsm_ln )

AUTO_TPACK( mat_trsm_ln )
{
	ADD_BLAS3_CASES_P3( trsm_ln, float )
	ADD_BLAS3_CASES_P3( trsm_ln, double )
}


DEF_BLAS3_CASES_P3_TR( trsm_lt )

AUTO_TPACK( mat_trsm_lt )
{
	ADD_BLAS3_CASES_P3( trsm_lt, float )
	ADD_BLAS3_CASES_P3( trsm_lt, double )
}

DEF_BLAS3_CASES_P3_TR( trsm_rn )

AUTO_TPACK( mat_trsm_rn )
{
	ADD_BLAS3_CASES_P3( trsm_rn, float )
	ADD_BLAS3_CASES_P3( trsm_rn, double )
}

DEF_BLAS3_CASES_P3_TR( trsm_rt )

AUTO_TPACK( mat_trsm_rt )
{
	ADD_BLAS3_CASES_P3( trsm_rt, float )
	ADD_BLAS3_CASES_P3( trsm_rt, double )
}

