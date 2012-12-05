/**
 * @file test_dense_acc.cpp
 *
 * Unit testing of accessors for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matexpr/dense_accessors.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary facilities

inline void fill_lin(dblock<double>& a)
{
	const index_t n = a.nelems();
	for (index_t i = 0; i < n; ++i)
	{
		a[i] = double(i+1);
	}
}


template<
	typename Tag1, typename Tag2,
	typename Acc, typename Ker, int M, int N>
void test_acc_eval()
{
	typedef typename mat_host<Tag1, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<Tag2, double, M, N>::mat_t dmat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_host<Tag1, double, M, N> src(m, n);
	mat_host<Tag2, double, M, N> dst(m, n);

	src.fill_lin();

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	typedef macc_policy<Acc, Ker> policy_t;
	typedef macc_scheme<policy_t, M, N> scheme_t;

	fill(dmat, 0.0);
	scheme_t sch = scheme_t::get_default(smat, dmat);
	sch.evaluate(smat, dmat);

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


#define TEST_ACC_EVAL( pname, acc, ker, stag, dtag ) \
	MN_CASE( pname, stag##_to_##dtag ) { \
		test_acc_eval<stag, dtag, acc, ker, M, N>(); } \
	BEGIN_TPACK( pname##_##stag##_to_##dtag ) \
		ADD_MN_CASE_3X3( pname, stag##_to_##dtag, DM, DN ) \
	END_TPACK

#define TEST_ACC_EVAL_V( pname, acc, ker, stag, dtag ) \
	MN_CASE( pname, stag##_to_##dtag ) { \
		test_acc_eval<stag, dtag, acc, ker, M, N>(); } \
	BEGIN_TPACK( pname##_##stag##_to_##dtag ) \
		ADD_MN_CASE( pname, stag##_to_##dtag, 1, 1 ) \
		ADD_MN_CASE( pname, stag##_to_##dtag, 0, 1 ) \
		ADD_MN_CASE( pname, stag##_to_##dtag, DM, 1 ) \
		ADD_MN_CASE( pname, stag##_to_##dtag, 1, 0 ) \
		ADD_MN_CASE( pname, stag##_to_##dtag, 1, DN ) \
	END_TPACK

#define TEST_LINEAR_SCALAR( stag, dtag ) \
	TEST_ACC_EVAL( linear_scalar, linear_macc, scalar_ker, stag, dtag )

#define TEST_LINEAR_SCALAR_V( stag, dtag ) \
	TEST_ACC_EVAL_V( linear_scalar, linear_macc, scalar_ker, stag, dtag )

#define ADD_LINEAR_SCALAR_PACK( stag, dtag ) \
	ADD_TPACK( linear_scalar_##stag##_to_##dtag )

#define TEST_PERCOL_SCALAR( stag, dtag ) \
	TEST_ACC_EVAL( percol_scalar, percol_macc, scalar_ker, stag, dtag )

#define ADD_PERCOL_SCALAR_PACK( stag, dtag ) \
	ADD_TPACK( percol_scalar_##stag##_to_##dtag )

TEST_LINEAR_SCALAR( cont, cont )
TEST_LINEAR_SCALAR( bloc, cont )
TEST_LINEAR_SCALAR( grid, cont )

TEST_LINEAR_SCALAR_V( cont, bloc )
TEST_LINEAR_SCALAR_V( bloc, bloc )
TEST_LINEAR_SCALAR_V( grid, bloc )

TEST_LINEAR_SCALAR_V( cont, grid )
TEST_LINEAR_SCALAR_V( bloc, grid )
TEST_LINEAR_SCALAR_V( grid, grid )

TEST_PERCOL_SCALAR( cont, cont )
TEST_PERCOL_SCALAR( bloc, cont )
TEST_PERCOL_SCALAR( grid, cont )

TEST_PERCOL_SCALAR( cont, bloc )
TEST_PERCOL_SCALAR( bloc, bloc )
TEST_PERCOL_SCALAR( grid, bloc )

TEST_PERCOL_SCALAR( cont, grid )
TEST_PERCOL_SCALAR( bloc, grid )
TEST_PERCOL_SCALAR( grid, grid )


BEGIN_MAIN_SUITE
	ADD_LINEAR_SCALAR_PACK( cont, cont )
	ADD_LINEAR_SCALAR_PACK( cont, bloc )
	ADD_LINEAR_SCALAR_PACK( cont, grid )
	ADD_LINEAR_SCALAR_PACK( bloc, cont )
	ADD_LINEAR_SCALAR_PACK( bloc, bloc )
	ADD_LINEAR_SCALAR_PACK( bloc, grid )
	ADD_LINEAR_SCALAR_PACK( grid, cont )
	ADD_LINEAR_SCALAR_PACK( grid, bloc )
	ADD_LINEAR_SCALAR_PACK( grid, grid )

	ADD_PERCOL_SCALAR_PACK( cont, cont )
	ADD_PERCOL_SCALAR_PACK( cont, bloc )
	ADD_PERCOL_SCALAR_PACK( cont, grid )
	ADD_PERCOL_SCALAR_PACK( bloc, cont )
	ADD_PERCOL_SCALAR_PACK( bloc, bloc )
	ADD_PERCOL_SCALAR_PACK( bloc, grid )
	ADD_PERCOL_SCALAR_PACK( grid, cont )
	ADD_PERCOL_SCALAR_PACK( grid, bloc )
	ADD_PERCOL_SCALAR_PACK( grid, grid )
END_MAIN_SUITE




