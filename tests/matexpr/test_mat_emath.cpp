/**
 * @file test_mat_emath.cpp
 *
 * Unit test of elementary functions on matrices
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_emath.h>

#include <cstdlib>
#include <functional>

using namespace lmat;
using namespace lmat::test;

const int DM = 8;
const int DN = 6;
const index_t LDim = 12;

template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}

template<class A, class B, typename T>
bool my_is_approx(const A& a, const B& b, const T& tol)
{
	if ( have_same_shape(a, b) )
	{
		index_t m = a.nrows();
		index_t n = a.ncolumns();

		return ltest::test_matrix_approx(m, n, a, b, tol);
	}
	else
	{
		return false;
	}
}

template<typename FTag, typename... Args>
inline void check_policy(const map_expr<FTag, Args...>& )
{
	typedef map_expr<FTag, Args...> Expr;
	typedef preferred_macc_policy<Expr> pmap;

	typedef typename matrix_traits<Expr>::value_type T;
	const int pw = (int)math::simd_traits<T, default_simd_kind>::pack_width;
	ASSERT_TRUE( pmap::prefer_linear );

	int M = meta::nrows<Expr>::value;
	int N = meta::ncols<Expr>::value;

	bool use_simd = ((M * N) % pw == 0) &&
			meta::has_simd_support<FTag, T, default_simd_kind>::value;
	ASSERT_EQ( pmap::prefer_simd, use_simd );
}

template<typename T, int M, int N>
struct tspec{ };

#define TEST_UNARY_EMATH(matfun, scafun, al, ah, tol) \
	MN_CASE( mat_emath, matfun ) { \
	typedef dense_matrix<double, M, N> mat_t; \
	const index_t m = M == 0 ? DM : M; \
	const index_t n = N == 0 ? DN : N; \
	mat_t A(m, n); fill_ran(A, al, ah); \
	check_policy( matfun(A) ); \
	mat_t R_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R_r[i] = scafun(A[i]); \
	mat_t R = matfun(A); \
	ASSERT_TRUE( my_is_approx(R, R_r, tol) ); \
	} \
	BEGIN_TPACK( mat_##matfun ) \
	ADD_MN_CASE_3X3( mat_emath, matfun, DM, DN ) \
	END_TPACK

#define TEST_BINARY_EMATH(matfun, scafun, al, ah, bl, bh, tol) \
	MN_CASE( mat_emath, matfun ) { \
	typedef dense_matrix<double, M, N> mat_t; \
	const index_t m = M == 0 ? DM : M; \
	const index_t n = N == 0 ? DN : N; \
	mat_t A(m, n); fill_ran(A, al, ah); \
	mat_t B(m, n); fill_ran(B, bl, bh); \
	check_policy( matfun(A, B) ); \
	mat_t R1_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R1_r[i] = scafun(A[i], B[i]); \
	mat_t R1 = matfun(A, B); \
	ASSERT_TRUE( my_is_approx(R1, R1_r, tol) ); \
	mat_t R2_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R2_r[i] = scafun(A[i], B[0]); \
	mat_t R2 = matfun(A, B[0]); \
	ASSERT_TRUE( my_is_approx(R2, R2_r, tol) ); \
	mat_t R3_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R3_r[i] = scafun(A[0], B[i]); \
	mat_t R3 = matfun(A[0], B); \
	ASSERT_TRUE( my_is_approx(R3, R3_r, tol) ); \
	} \
	BEGIN_TPACK( mat_##matfun ) \
	ADD_MN_CASE_3X3( mat_emath, matfun, DM, DN ) \
	END_TPACK

#define TEST_TERNARY_EMATH(matfun, scafun, al, ah, bl, bh, cl, ch, tol) \
	MN_CASE( mat_emath, matfun ) { \
	typedef dense_matrix<double, M, N> mat_t; \
	const index_t m = M == 0 ? DM : M; \
	const index_t n = N == 0 ? DN : N; \
	mat_t A(m, n); fill_ran(A, al, ah); \
	mat_t B(m, n); fill_ran(B, bl, bh); \
	mat_t C(m, n); fill_ran(C, cl, ch); \
	check_policy( matfun(A, B, C) ); \
	mat_t R1_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R1_r[i] = scafun(A[i], B[i], C[i]); \
	mat_t R1 = matfun(A, B, C); \
	ASSERT_TRUE( my_is_approx(R1, R1_r, tol) ); \
	mat_t R2_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R2_r[i] = scafun(A[0], B[i], C[i]); \
	mat_t R2 = matfun(A[0], B, C); \
	ASSERT_TRUE( my_is_approx(R2, R2_r, tol) ); \
	mat_t R3_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R3_r[i] = scafun(A[i], B[0], C[i]); \
	mat_t R3 = matfun(A, B[0], C); \
	ASSERT_TRUE( my_is_approx(R3, R3_r, tol) ); \
	mat_t R4_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R4_r[i] = scafun(A[i], B[i], C[0]); \
	mat_t R4 = matfun(A, B, C[0]); \
	ASSERT_TRUE( my_is_approx(R4, R4_r, tol) ); \
	mat_t R5_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R5_r[i] = scafun(A[0], B[0], C[i]); \
	mat_t R5 = matfun(A[0], B[0], C); \
	ASSERT_TRUE( my_is_approx(R5, R5_r, tol) ); \
	mat_t R6_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R6_r[i] = scafun(A[0], B[i], C[0]); \
	mat_t R6 = matfun(A[0], B, C[0]); \
	ASSERT_TRUE( my_is_approx(R6, R6_r, tol) ); \
	mat_t R7_r(m, n); \
	for (index_t i = 0; i < m * n; ++i) R7_r[i] = scafun(A[i], B[0], C[0]); \
	mat_t R7 = matfun(A, B[0], C[0]); \
	ASSERT_TRUE( my_is_approx(R7, R7_r, tol) ); \
	} \
	BEGIN_TPACK( mat_##matfun ) \
	ADD_MN_CASE_3X3( mat_emath, matfun, DM, DN ) \
	END_TPACK


template<typename T>
inline T my_sqr(T x) { return x * x; }

template<typename T>
inline T my_cube(T x) { return x * x * x; }

template<typename T>
inline T my_rcp(T x) { return T(1) / x; }

template<typename T>
inline T my_rsqrt(T x) { return T(1) / std::sqrt(x); }

template<typename T>
inline T my_fma(T x, T y, T z) { return x * y + z; }

template<typename T>
inline T my_clamp(T x, T y, T z)
{
	T r = x;
	if (r < y) r = y;
	if (r > z) r = z;
	return r;
}


TEST_TERNARY_EMATH(fma, my_fma, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, 1.0e-15)
TEST_TERNARY_EMATH(clamp, my_clamp, -1.0, 1.0, -1.0, 0.0, 0.0, 1.0, 1.0e-15)

TEST_UNARY_EMATH(sqrt, std::sqrt, 0.0, 5.0, 1.0e-14)
TEST_UNARY_EMATH(rcp, my_rcp, 0.5, 5.0, 1.0e-14)
TEST_UNARY_EMATH(rsqrt, my_rsqrt, 0.5, 5.0, 1.0e-14)
TEST_BINARY_EMATH(pow, std::pow, 0.5, 2.0, 0.5, 2.0, 1.0e-14)

TEST_UNARY_EMATH(floor, std::floor, -5.0, 5.0, 1.0e-16)
TEST_UNARY_EMATH(ceil, std::ceil, -5.0, 5.0, 1.0e-16)

TEST_UNARY_EMATH(exp, std::exp, -1.0, 3.0, 1.0e-14)
TEST_UNARY_EMATH(log, std::log, 0.5, 10.0, 1.0e-14)
TEST_UNARY_EMATH(log10, std::log10, 0.5, 20.0, 1.0e-14)
TEST_UNARY_EMATH(xlogx, math::xlogx, -1.0, 3.0, 1.0e-14)
TEST_BINARY_EMATH(xlogy, math::xlogy, -1.0, 3.0, 1.0, 4.0, 1.0e-14)

TEST_UNARY_EMATH(sin, std::sin, -3.0, 3.0, 1.0e-14)
TEST_UNARY_EMATH(cos, std::cos, -3.0, 3.0, 1.0e-14)
TEST_UNARY_EMATH(tan, std::tan, -1.0, 1.0, 1.0e-14)

TEST_UNARY_EMATH(asin, std::asin, -1.0, 1.0, 1.0e-14)
TEST_UNARY_EMATH(acos, std::acos, -1.0, 1.0, 1.0e-14)
TEST_UNARY_EMATH(atan, std::atan, -1.0, 1.0, 1.0e-14)
TEST_BINARY_EMATH(atan2, std::atan2, -3.0, 3.0, -3.0, 3.0, 1.0e-14)

TEST_UNARY_EMATH(sinh, std::sinh, -2.0, 2.0, 1.0e-14)
TEST_UNARY_EMATH(cosh, std::cosh, -2.0, 2.0, 1.0e-14)
TEST_UNARY_EMATH(tanh, std::tanh, -2.0, 2.0, 1.0e-14)

#ifdef LMAT_HAS_CXX11_MATH

TEST_UNARY_EMATH(cbrt, std::cbrt, -10.0, 10.0, 1.0e-14)
TEST_BINARY_EMATH(hypot, std::hypot, -3.0, 3.0, -3.0, 3.0, 1.0e-14)

TEST_UNARY_EMATH(round, std::round, -5.0, 5.0, 1.0e-16)
TEST_UNARY_EMATH(trunc, std::trunc, -5.0, 5.0, 1.0e-16)

TEST_UNARY_EMATH(exp2, std::exp2, -1.0, 3.0, 1.0e-14)
TEST_UNARY_EMATH(log2, std::log2, 0.5, 10.0, 1.0e-14)
TEST_UNARY_EMATH(expm1, std::expm1, -1.0, 1.0, 1.0e-15)
TEST_UNARY_EMATH(log1p, std::log1p, 0.1, 1.0, 1.0e-15)

TEST_UNARY_EMATH(asinh, std::asinh, -2.0, 2.0, 1.0e-14)
TEST_UNARY_EMATH(acosh, std::acosh, 1.0, 5.0, 1.0e-14)
TEST_UNARY_EMATH(atanh, std::atanh, -1.0, 1.0, 1.0e-14)

TEST_UNARY_EMATH(erf, std::erf, -2.0, 2.0, 1.0e-14)
TEST_UNARY_EMATH(erfc, std::erfc, -2.0, 2.0, 1.0e-14)
TEST_UNARY_EMATH(lgamma, std::lgamma, 1.0, 5.0, 1.0e-12)
TEST_UNARY_EMATH(tgamma, std::tgamma, 1.0, 3.0, 1.0e-12)

#endif

TEST_UNARY_EMATH(norminv, math::norminv, 0.0, 1.0, 1.0e-12)


BEGIN_MAIN_SUITE

	ADD_TPACK( mat_fma )
	ADD_TPACK( mat_clamp )

	ADD_TPACK( mat_sqrt )
	ADD_TPACK( mat_rcp )
	ADD_TPACK( mat_rsqrt )
	ADD_TPACK( mat_pow )

	ADD_TPACK( mat_floor )
	ADD_TPACK( mat_ceil )

	ADD_TPACK( mat_exp )
	ADD_TPACK( mat_log )
	ADD_TPACK( mat_log10 )
	ADD_TPACK( mat_xlogy )
	ADD_TPACK( mat_xlogx )

	ADD_TPACK( mat_sin )
	ADD_TPACK( mat_cos )
	ADD_TPACK( mat_tan )

	ADD_TPACK( mat_asin )
	ADD_TPACK( mat_acos )
	ADD_TPACK( mat_atan )
	ADD_TPACK( mat_atan2 )

	ADD_TPACK( mat_sinh )
	ADD_TPACK( mat_cosh )
	ADD_TPACK( mat_tanh )

#ifdef LMAT_HAS_CXX11_MATH

	ADD_TPACK( mat_cbrt )
	ADD_TPACK( mat_hypot )

	ADD_TPACK( mat_round )
	ADD_TPACK( mat_trunc )

	ADD_TPACK( mat_exp2 )
	ADD_TPACK( mat_log2 )
	ADD_TPACK( mat_expm1 )
	ADD_TPACK( mat_log1p )

	ADD_TPACK( mat_asinh )
	ADD_TPACK( mat_acosh )
	ADD_TPACK( mat_atanh )

	ADD_TPACK( mat_erf )
	ADD_TPACK( mat_erfc )
	ADD_TPACK( mat_lgamma )
	ADD_TPACK( mat_tgamma )
#endif

	ADD_TPACK( mat_norminv )
END_MAIN_SUITE





