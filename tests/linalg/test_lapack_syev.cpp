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
using lmat::lapack::syevr;

using lmat::lapack::evr_I;
using lmat::lapack::evr_V;


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


template<typename W>
void regularize_sign(W& a)
{
	index_t m = a.nrows();
	index_t n = a.ncolumns();

	for (index_t j = 0; j < n; ++j)
	{
		if (a(0, j) < 0)
		{
			for (index_t i = 0; i < m; ++i)
				a(i, j) = - a(i, j);
		}
	}
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


T_CASE( mat_syev )
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

	T tol0 = (T)(sizeof(T) == 4 ? 1.0e-4 : 1.0e-13);
	T tol = (T)(sizeof(T) == 4 ? 1.0e-4 : 1.0e-10);

	test_syev(a, w0, w, V, tol0, tol);
}


T_CASE( mat_syevd )
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


T_CASE( mat_syevr )
{
	typedef mat_host<bloc, T, 0, 0> host_t;
	typedef typename host_t::mat_t mat_t;
	typedef typename host_t::cmat_t cmat_t;

	index_t m = DM;
	host_t a_host(m, m);
	mat_t a = a_host.get_mat();
	fill_rand_pdm(a);

	index_t ret;

	dense_col<T> w0;
	ret = syevr(a, w0);

	ASSERT_EQ( ret, m );
	ASSERT_EQ( w0.nrows(), m );
	ASSERT_EQ( w0.ncolumns(), 1 );

	dense_col<T> w;
	dense_matrix<T> V;
	ret = syevr(a, w, V);

	ASSERT_EQ( ret, m );
	ASSERT_EQ( w.nrows(), m );
	ASSERT_EQ( w.ncolumns(), 1 );
	ASSERT_EQ( V.nrows(), m );
	ASSERT_EQ( V.ncolumns(), m );

	T tol0 = (T)(sizeof(T) == 4 ? 1.0e-5 : 1.0e-12);
	T tol = (T)(sizeof(T) == 4 ? 1.0e-4 : 1.0e-10);

	regularize_sign(V);

	test_syev(a, w0, w, V, tol0, tol);

	// select index range

	dense_col<T> wI;
	ret = syevr(a, wI, evr_I(1, m-1));

	ASSERT_EQ( ret, m-2 );
	ASSERT_EQ( wI.nrows(), m-2 );
	ASSERT_EQ( wI.ncolumns(), 1 );

	auto w0I = w0(range(1, m-2));
	ASSERT_VEC_APPROX( m-2, wI, w0I, tol0 );

	dense_matrix<T> VI;
	ret = syevr(a, wI, VI, evr_I(1, m-1));

	ASSERT_EQ( ret, m-2 );
	ASSERT_EQ( wI.nrows(), m-2 );
	ASSERT_EQ( wI.ncolumns(), 1 );
	ASSERT_EQ( VI.nrows(), m );
	ASSERT_EQ( VI.ncolumns(), m-2 );

	auto V0I = V(whole(), range(1, m-2));
	regularize_sign(VI);

	ASSERT_VEC_APPROX( m-2, wI, w0I, tol0 );
	ASSERT_MAT_APPROX( m, m-2, VI, V0I, tol0 );

	// select value range

	T vl = (w[0] + w[1]) * T(0.5);
	T vu = (w[m-2] + w[m-1]) * T(0.5);

	dense_col<T> wV;
	ret = syevr(a, wV, evr_V(vl, vu));

	ASSERT_EQ( ret, m-2 );
	ASSERT_EQ( wV.nrows(), m );
	ASSERT_EQ( wV.ncolumns(), 1 );

	ASSERT_VEC_APPROX( m-2, wV, w0I, tol );

	dense_matrix<T> VV;
	ret = syevr(a, wV, VV, evr_V(vl, vu));

	ASSERT_EQ( ret, m-2 );
	ASSERT_EQ( wV.nrows(), m );
	ASSERT_EQ( wV.ncolumns(), 1 );
	ASSERT_EQ( VV.nrows(), m );
	ASSERT_EQ( VV.ncolumns(), m );

	regularize_sign(VV);

	ASSERT_VEC_APPROX( m-2, wV, w0I, tol0 );
	ASSERT_MAT_APPROX( m, m-2, VV, V0I, tol0 );
}


AUTO_TPACK( mat_syev )
{
	ADD_T_CASE( mat_syev, float )
	ADD_T_CASE( mat_syev, double )
}

AUTO_TPACK( mat_syevd )
{
	ADD_T_CASE( mat_syevd, float )
	ADD_T_CASE( mat_syevd, double )
}

AUTO_TPACK( mat_syevr )
{
	ADD_T_CASE( mat_syevr, float )
	ADD_T_CASE( mat_syevr, double )
}



