/**
 * @file test_lapack_svd.cpp
 *
 * @brief Unit testing of Singular Value Decomposition
 *
 * @author Dahua Lin
 */


#include "linalg_test_base.h"
#include <light_mat/linalg/blas_l3.h>
#include <light_mat/linalg/lapack_svd.h>

using namespace lmat;
using namespace lmat::test;

using lmat::lapack::gesvd;
using lmat::lapack::gesdd;


#define DISP(a) std::printf("\n" #a "=\n"); printf_mat("%10.4g ", a); std::printf("\n")

template<typename T, class Mat>
bool is_orth( const IRegularMatrix<Mat, T>& a, char trans, T tol)
{
	index_t m = a.nrows();
	index_t n = a.ncolumns();

	if (trans == 'N' || trans == 'n')
	{
		dense_matrix<T> e(n, n, zero());
		for (index_t i = 0; i < n; ++i) e(i, i) = T(1);

		dense_matrix<T> r(n, n, zero());
		blas::gemm(a, a, r, 'T', 'N');

		return ltest::test_matrix_approx(n, n, r, e, tol);
	}
	else // 'T'
	{
		dense_matrix<T> e(m, m, zero());
		for (index_t i = 0; i < m; ++i) e(i, i) = T(1);

		dense_matrix<T> r(m, m, zero());
		blas::gemm(a, a, r, 'N', 'T');

		return ltest::test_matrix_approx(m, m, r, e, tol);
	}
}


template<typename A>
bool is_sorted_desc(const A& a)
{
	index_t n = a.nelems();
	for (index_t i = 0; i < n-1; ++i)
	{
		if (a[i] < a[i+1]) return false;
	}

	return true;
}


template<typename T, typename A, typename U, typename S, typename VT>
bool check_svd(const IRegularMatrix<A, T>& a, const IRegularMatrix<S, T>& s,
		const IRegularMatrix<U, T>& u, const IRegularMatrix<VT, T>& vt, T tol)
{
	index_t m = a.nrows();
	index_t n = a.ncolumns();

	index_t sm = u.ncolumns();
	index_t sn = vt.nrows();

	const S& s_ = s.derived();

	dense_matrix<T> smat(sm, sn, zero());
	for (index_t i = 0; i < math::min(sm, sn); ++i)
	{
		smat(i, i) = s_[i];
	}

	dense_matrix<T> us(m, sn);
	blas::gemm(u, smat, us);

	dense_matrix<T> prod(m, n);
	blas::gemm(us, vt, prod);

	return ltest::test_matrix_approx(m, n, a, prod, tol);
}


template<typename T>
void test_svd( index_t m, index_t n )
{
	dense_matrix<T> a(m, n);
	do_fill_rand(a.ptr_data(), m * n);

	index_t k = math::min(m, n);
	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);

	dense_col<T> s0;

	gesvd(a, s0);
	ASSERT_EQ(s0.nrows(), k);
	ASSERT_EQ(s0.ncolumns(), 1);
	ASSERT_TRUE( is_sorted_desc(s0) );

	dense_col<T> s;
	dense_matrix<T> u;
	dense_matrix<T> vt;

	gesvd(a, s, u, vt);

	ASSERT_EQ(s.nrows(), k);
	ASSERT_EQ(s.ncolumns(), 1);
	ASSERT_EQ(u.nrows(), m);
	ASSERT_EQ(u.ncolumns(), m);
	ASSERT_EQ(vt.nrows(), n);
	ASSERT_EQ(vt.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s) );
	ASSERT_VEC_APPROX(k, s, s0, tol);

	ASSERT_TRUE( is_orth(u, 'N', tol) );
	ASSERT_TRUE( is_orth(vt, 'T', tol) );

	ASSERT_TRUE( check_svd(a, s, u, vt, tol) );

	dense_col<T> s1;
	dense_matrix<T> u1;
	dense_matrix<T> vt1;

	gesvd(a, s1, u1, vt1, 'A', 'S');

	ASSERT_EQ(s1.nrows(), k);
	ASSERT_EQ(s1.ncolumns(), 1);
	ASSERT_EQ(u1.nrows(), m);
	ASSERT_EQ(u1.ncolumns(), m);
	ASSERT_EQ(vt1.nrows(), k);
	ASSERT_EQ(vt1.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s1) );
	ASSERT_VEC_APPROX(k, s1, s0, tol);

	ASSERT_TRUE( is_orth(u1, 'N', tol) );
	ASSERT_TRUE( is_orth(vt1, 'T', tol) );
	ASSERT_TRUE( check_svd(a, s1, u1, vt1, tol) );

	dense_col<T> s2;
	dense_matrix<T> u2;
	dense_matrix<T> vt2;

	gesvd(a, s2, u2, vt2, 'S', 'A');

	ASSERT_EQ(s2.nrows(), k);
	ASSERT_EQ(s2.ncolumns(), 1);
	ASSERT_EQ(u2.nrows(), m);
	ASSERT_EQ(u2.ncolumns(), k);
	ASSERT_EQ(vt2.nrows(), n);
	ASSERT_EQ(vt2.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s2) );
	ASSERT_VEC_APPROX(k, s2, s0, tol);
	ASSERT_TRUE( is_orth(u2, 'N', tol) );
	ASSERT_TRUE( is_orth(vt2, 'T', tol) );

	ASSERT_TRUE( check_svd(a, s2, u2, vt2, tol) );

	dense_col<T> s3;
	dense_matrix<T> u3;
	dense_matrix<T> vt3;

	gesvd(a, s3, u3, vt3, 'S', 'S');

	ASSERT_EQ(s3.nrows(), k);
	ASSERT_EQ(s3.ncolumns(), 1);
	ASSERT_EQ(u3.nrows(), m);
	ASSERT_EQ(u3.ncolumns(), k);
	ASSERT_EQ(vt3.nrows(), k);
	ASSERT_EQ(vt3.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s3) );
	ASSERT_VEC_APPROX(k, s3, s0, tol);

	ASSERT_TRUE( is_orth(u3, 'N', tol) );
	ASSERT_TRUE( is_orth(vt3, 'T', tol) );

	ASSERT_TRUE( check_svd(a, s3, u3, vt3, tol) );

	dense_col<T> s4;
	dense_matrix<T> u4;
	dense_matrix<T> vt4;

	gesvd(a, s4, u4, vt4, 'A', 'N');
	ASSERT_EQ( u4.nrows(), m );
	ASSERT_EQ( u4.ncolumns(), m );
	ASSERT_MAT_APPROX(m, m, u4, u, tol);

	dense_col<T> s5;
	dense_matrix<T> u5;
	dense_matrix<T> vt5;

	gesvd(a, s5, u5, vt5, 'S', 'N');
	ASSERT_EQ( u5.nrows(), m );
	ASSERT_EQ( u5.ncolumns(), k );
	ASSERT_MAT_APPROX(m, k, u5, u2, tol);

	dense_col<T> s6;
	dense_matrix<T> u6;
	dense_matrix<T> vt6;

	gesvd(a, s6, u6, vt6, 'N', 'A');
	ASSERT_EQ( vt6.nrows(), n );
	ASSERT_EQ( vt6.ncolumns(), n );
	ASSERT_MAT_APPROX(n, n, vt6, vt, tol);

	dense_col<T> s7;
	dense_matrix<T> u7;
	dense_matrix<T> vt7;

	gesvd(a, s7, u7, vt7, 'N', 'S');
	ASSERT_EQ( vt7.nrows(), k );
	ASSERT_EQ( vt7.ncolumns(), n );
	ASSERT_MAT_APPROX(k, n, vt7, vt1, tol);
}


template<typename T>
void test_sdd( index_t m, index_t n )
{
	dense_matrix<T> a(m, n);
	do_fill_rand(a.ptr_data(), m * n);

	index_t k = math::min(m, n);
	T tol = (T)(sizeof(T) == 4 ? 2.0e-5 : 1.0e-10);

	dense_col<T> s0;

	gesvd(a, s0);
	ASSERT_EQ(s0.nrows(), k);
	ASSERT_EQ(s0.ncolumns(), 1);
	ASSERT_TRUE( is_sorted_desc(s0) );

	dense_col<T> s;
	dense_matrix<T> u;
	dense_matrix<T> vt;

	gesdd(a, s, u, vt);

	ASSERT_EQ(s.nrows(), k);
	ASSERT_EQ(s.ncolumns(), 1);
	ASSERT_EQ(u.nrows(), m);
	ASSERT_EQ(u.ncolumns(), m);
	ASSERT_EQ(vt.nrows(), n);
	ASSERT_EQ(vt.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s) );
	ASSERT_VEC_APPROX(k, s, s0, tol);

	ASSERT_TRUE( is_orth(u, 'N', tol) );
	ASSERT_TRUE( is_orth(vt, 'T', tol) );

	ASSERT_TRUE( check_svd(a, s, u, vt, tol) );

	dense_col<T> s3;
	dense_matrix<T> u3;
	dense_matrix<T> vt3;

	gesdd(a, s3, u3, vt3, 'S');

	ASSERT_EQ(s3.nrows(), k);
	ASSERT_EQ(s3.ncolumns(), 1);
	ASSERT_EQ(u3.nrows(), m);
	ASSERT_EQ(u3.ncolumns(), k);
	ASSERT_EQ(vt3.nrows(), k);
	ASSERT_EQ(vt3.ncolumns(), n);

	ASSERT_TRUE( is_sorted_desc(s3) );
	ASSERT_VEC_APPROX(k, s3, s0, tol);

	ASSERT_TRUE( is_orth(u3, 'N', tol) );
	ASSERT_TRUE( is_orth(vt3, 'T', tol) );

	ASSERT_TRUE( check_svd(a, s3, u3, vt3, tol) );
}



T_CASE( mat_svd, svd_eq )
{
	test_svd<T>(8, 8);
}

T_CASE( mat_svd, svd_gt )
{
	test_svd<T>(8, 5);
}

T_CASE( mat_svd, svd_lt )
{
	test_svd<T>(5, 8);
}


T_CASE( mat_svd, sdd_eq )
{
	test_sdd<T>(8, 8);
}

T_CASE( mat_svd, sdd_gt )
{
	test_sdd<T>(8, 5);
}

T_CASE( mat_svd, sdd_lt )
{
	test_sdd<T>(5, 8);
}


BEGIN_TPACK( mat_svd )
	ADD_T_CASE( mat_svd, svd_eq, float )
	ADD_T_CASE( mat_svd, svd_gt, float )
	ADD_T_CASE( mat_svd, svd_lt, float )
	ADD_T_CASE( mat_svd, svd_eq, double )
	ADD_T_CASE( mat_svd, svd_gt, double )
	ADD_T_CASE( mat_svd, svd_lt, double )
END_TPACK

BEGIN_TPACK( mat_sdd )
	ADD_T_CASE( mat_svd, sdd_eq, float )
	ADD_T_CASE( mat_svd, sdd_gt, float )
	ADD_T_CASE( mat_svd, sdd_lt, float )
	ADD_T_CASE( mat_svd, sdd_eq, double )
	ADD_T_CASE( mat_svd, sdd_gt, double )
	ADD_T_CASE( mat_svd, sdd_lt, double )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_svd )
	ADD_TPACK( mat_sdd )
END_MAIN_SUITE



