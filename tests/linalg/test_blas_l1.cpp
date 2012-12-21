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


TMN_CASE( mat_blas1, asum_cont ) { test_asum<cont, T, M, N>(); }


BEGIN_TPACK( mat_asum_cont )
	ADD_TMN_CASE( mat_blas1, asum_cont, float, 0, 0 )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_asum_cont )
END_MAIN_SUITE

