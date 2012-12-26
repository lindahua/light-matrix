/**
 * @file test_mat_pred.cpp
 *
 * @brief Unit testing of element-wise predicates on matrices
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_pred.h>

using namespace lmat;
using namespace lmat::test;

typedef mask_t<double> fmsk;

const int DM = 9;
const int DN = 8;


template<class A, class B>
bool my_is_equal(const A& a, const B& b)
{
	if ( have_same_shape(a, b) )
	{
		index_t m = a.nrows();
		index_t n = a.ncolumns();

		return ltest::test_matrix_equal(m, n, a, b);
	}
	else
	{
		return false;
	}
}


#define DEF_NUMCOMP_CASE( name, op, aexpr, bexpr, cexpr ) \
		MN_CASE( mat_comp, name ) { \
			typedef dense_matrix<double, M, N> mat_t; \
			typedef dense_matrix<bool, M, N> bmat_t; \
			typedef dense_matrix<fmsk, M, N> mmat_t; \
			const index_t m = M == 0 ? DM : M; \
			const index_t n = N == 0 ? DN : N; \
			mat_t A(m, n); \
			mat_t B(m, n); \
			double c = double( cexpr ); \
			for (index_t i = 0; i < m * n; ++i) A[i] = double( aexpr ); \
			for (index_t i = 0; i < m * n; ++i) B[i] = double( bexpr ); \
			bmat_t AB_r(m, n); \
			for (index_t i = 0; i < m * n; ++i) AB_r[i] = (A[i] op B[i]); \
			bmat_t AC_r(m, n); \
			for (index_t i = 0; i < m * n; ++i) AC_r[i] = (A[i] op c); \
			bmat_t CB_r(m, n); \
			for (index_t i = 0; i < m * n; ++i) CB_r[i] = (c op B[i]); \
			mmat_t ABm = (A op B); \
			bmat_t AB = to_bool(ABm); \
			ASSERT_TRUE( my_is_equal(AB, AB_r) ); \
			mmat_t ACm = (A op c); \
			bmat_t AC = to_bool(ACm); \
			ASSERT_TRUE( my_is_equal(AC, AC_r) ); \
			mmat_t CBm = (c op B); \
			bmat_t CB = to_bool(CBm); \
			ASSERT_TRUE( my_is_equal(CB, CB_r) ); }

#define DEF_NUMCOMP_TESTS( name, op ) \
		DEF_NUMCOMP_CASE( name, op, i+1, (i+1) * (i%2==1), (m*n)/2 ) \
		BEGIN_TPACK( mat_##name ) \
			ADD_MN_CASE_3X3( mat_comp, name, DM, DN ) \
		END_TPACK


DEF_NUMCOMP_TESTS( eq, == )
DEF_NUMCOMP_TESTS( ne, != )
DEF_NUMCOMP_TESTS( ge, >= )
DEF_NUMCOMP_TESTS( gt, >  )
DEF_NUMCOMP_TESTS( le, <= )
DEF_NUMCOMP_TESTS( lt, <  )

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_eq )
	ADD_TPACK( mat_ne )
	ADD_TPACK( mat_ge )
	ADD_TPACK( mat_gt )
	ADD_TPACK( mat_le )
	ADD_TPACK( mat_lt )
END_MAIN_SUITE



