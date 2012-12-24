/**
 * @file test_lapack_syev.cpp
 *
 * @brief Unit testing of symmetric eigenvalue problems
 *
 * @author Dahua Lin
 */

#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l3.h>
#include <light_mat/linalg/lapack_syev.h>

using namespace lmat;
using namespace lmat::test;

using lmat::lapack::syev;
using lmat::lapack::syevd;


template<typename W>
bool is_sorted_asc_(const W& w)
{
	index_t n = w.nelems();
	for (index_t i = 0; i < n-1; ++i)
	{
		if (w[i] > w[i+1]) return false;
	}

	return true;
}

template<typename T, class A, class CW, class CV>
void test_syev(const A& a, const CW& w0, const CW& w, const CV& V, T tol0, T tol)
{
	ASSERT_TRUE( is_sorted_asc_(w0) );
	ASSERT_TRUE( is_sorted_asc_(w) );

	index_t m = a.nrows();
	ASSERT_VEC_APPROX(m, w0, w, tol0);

	dense_matrix<T> D(m, m, zero());
	for (index_t i = 0; i < m; ++i) D(i, i) = w[i];

	dense_matrix<T> L(m, m, zero());
	dense_matrix<T> R(m, m, zero());

	blas::gemm(a, V, L);
	blas::gemm(V, D, R);
	ASSERT_MAT_APPROX(m, m, L, R, tol);
}


T_CASE( mat_syev, syev )
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::mat_t mat_t;
	typedef typename host_t::cmat_t cmat_t;

	index_t m = DM;
	host_t a_host(m, m);
	mat_t a = a_host.get_mat();
	fill_rand_pdm(a);

	dense_col<T> w0(m);
	syev(a, w0);

	dense_col<T> w(m);
	dense_matrix<T> V(m, m);

	syev(a, w, V);

	T tol0 = (T)(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	T tol = (T)(sizeof(T) == 4 ? 1.0e-4 : 1.0e-10);

	test_syev(a, w0, w, V, tol0, tol);
}


T_CASE( mat_syev, syevd )
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::mat_t mat_t;
	typedef typename host_t::cmat_t cmat_t;

	index_t m = DM;
	host_t a_host(m, m);
	mat_t a = a_host.get_mat();
	fill_rand_pdm(a);

	dense_col<T> w0(m);
	syevd(a, w0);

	dense_col<T> w(m);
	dense_matrix<T> V(m, m);

	syevd(a, w, V);

	T tol0 = (T)(sizeof(T) == 4 ? 1.0e-6 : 1.0e-13);
	T tol = (T)(sizeof(T) == 4 ? 1.0e-4 : 1.0e-10);

	test_syev(a, w0, w, V, tol0, tol);
}


BEGIN_TPACK( mat_syev )
	ADD_T_CASE( mat_syev, syev, float )
	ADD_T_CASE( mat_syev, syev, double )
END_TPACK

BEGIN_TPACK( mat_syevd )
	ADD_T_CASE( mat_syev, syevd, float )
	ADD_T_CASE( mat_syev, syevd, double )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_syev )
	ADD_TPACK( mat_syevd )
END_MAIN_SUITE


