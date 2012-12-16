/**
 * @file test_mat_arith.cpp
 *
 * Unit testing for Basic matrix arithmetics
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include "multimat_supp.h"

#include <light_mat/matexpr/matrix_arith.h>
#include <light_mat/matexpr/matrix_ewise_eval.h>

using namespace lmat;
using namespace lmat::test;

typedef dense_matrix<double> dmat_t;


template<int M, int N>
void fill_ran(dense_matrix<double, M, N>& X, double a, double b)
{
	for (index_t i = 0; i < X.nelems(); ++i)
	{
		X[i] = a + (double(std::rand()) / RAND_MAX) * (b - a);
	}
}


template<class Tag1, class Tag2, int M, int N>
void test_scheme_choice()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef typename mat_host<Tag1, double, M, N>::cmat_t cmat1_t;
	typedef typename mat_host<Tag2, double, M, N>::cmat_t cmat2_t;

	mat_host<Tag1, double, M, N> a_src(m, n);
	mat_host<Tag2, double, M, N> b_src(m, n);

	cmat1_t a = a_src.get_cmat();
	cmat2_t b = b_src.get_cmat();

	bool supp_lin1 = meta::is_continuous<cmat1_t>::value || meta::is_vector<cmat1_t>::value;
	bool supp_lin2 = meta::is_continuous<cmat2_t>::value || meta::is_vector<cmat2_t>::value;

	dmat_t r(m, n);

	bool expect_lin_u = get_default_macc_scheme(-a, r).use_linear();

	ASSERT_EQ( expect_lin_u, supp_lin1 );

	bool expect_lin_b = get_default_macc_scheme(a + b, r).use_linear();

	ASSERT_EQ( expect_lin_b, supp_lin1 && supp_lin2 );
}


MN_CASE( sch_choice, mat_mat )
{
	test_scheme_choice<cont, cont, M, N>();
}

MN_CASE( sch_choice, mat_blk )
{
	test_scheme_choice<cont, bloc, M, N>();
}

MN_CASE( sch_choice, mat_grid )
{
	test_scheme_choice<cont, grid, M, N>();
}

MN_CASE( sch_choice, blk_mat )
{
	test_scheme_choice<bloc, cont, M, N>();
}

MN_CASE( sch_choice, blk_blk )
{
	test_scheme_choice<bloc, bloc, M, N>();
}

MN_CASE( sch_choice, blk_grid )
{
	test_scheme_choice<bloc, grid, M, N>();
}

MN_CASE( sch_choice, grid_mat )
{
	test_scheme_choice<grid, cont, M, N>();
}

MN_CASE( sch_choice, grid_blk )
{
	test_scheme_choice<grid, bloc, M, N>();
}

MN_CASE( sch_choice, grid_grid )
{
	test_scheme_choice<grid, grid, M, N>();
}


MN_CASE( mat_arith, add )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c + B[i];

	// default evaluation

	mat_t AB = A + B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A + c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c + B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

}


MN_CASE( mat_arith, add_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] + B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] + c;

	// default evaluation

	mat_t AB(A);
	AB += B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC += c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}



MN_CASE( mat_arith, sub )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c - B[i];

	// default evaluation

	mat_t AB = A - B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A - c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c - B;
	ASSERT_TRUE( is_equal(CB, CB_r) );

}

MN_CASE( mat_arith, sub_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] - B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] - c;

	// default evaluation

	mat_t AB(A);
	AB -= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC -= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, mul )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c * B[i];

	// default evaluation

	mat_t AB = A * B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A * c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c * B;
	ASSERT_TRUE( is_equal(CB, CB_r) );
}

MN_CASE( mat_arith, mul_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 7.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(2 * i + 3);

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] * B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] * c;

	// default evaluation

	mat_t AB(A);
	AB *= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC *= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, div )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c / B[i];

	// default evaluation

	mat_t AB = A / B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = A / c;
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = c / B;
	ASSERT_TRUE( is_equal(CB, CB_r) );
}


MN_CASE( mat_arith, div_ip )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	double c = 4.0;

	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);
	for (index_t i = 0; i < m * n; ++i) B[i] = double(1 << (i % 5));

	// prepare ground-truth

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] / B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] / c;

	// default evaluation

	mat_t AB(A);
	AB /= B;
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC(A);
	AC /= c;
	ASSERT_TRUE( is_equal(AC, AC_r) );
}


MN_CASE( mat_arith, neg )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = -A[i];

	// default evaluation

	mat_t R = -A;
	ASSERT_TRUE( is_equal(R, R_r) );
}


MN_CASE( mat_arith, abs )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i+1) * (i % 3 - 1));

	// prepare ground-truth

	mat_t R_r(m, n);
	for (index_t i = 0; i < m * n; ++i) R_r[i] = std::abs(A[i]);

	// default evaluation

	mat_t R = abs(A);
	ASSERT_TRUE( is_equal(R, R_r) );
}

MN_CASE( mat_arith, max )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n); fill_ran(A, 0.0, 10.0);
	mat_t B(m, n); fill_ran(B, 0.0, 10.0);
	double c = 5.0;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] > B[i] ? A[i] : B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] > c ? A[i] : c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c > B[i] ? c : B[i];

	mat_t AB = (max)(A, B);
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = (max)(A, c);
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = (max)(c, B);
	ASSERT_TRUE( is_equal(CB, CB_r) );
}

MN_CASE( mat_arith, min )
{
	typedef dense_matrix<double, M, N> mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n); fill_ran(A, 0.0, 10.0);
	mat_t B(m, n); fill_ran(B, 0.0, 10.0);
	double c = 5.0;

	mat_t AB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AB_r[i] = A[i] < B[i] ? A[i] : B[i];

	mat_t AC_r(m, n);
	for (index_t i = 0; i < m * n; ++i) AC_r[i] = A[i] < c ? A[i] : c;

	mat_t CB_r(m, n);
	for (index_t i = 0; i < m * n; ++i) CB_r[i] = c < B[i] ? c : B[i];

	mat_t AB = (min)(A, B);
	ASSERT_TRUE( is_equal(AB, AB_r) );

	mat_t AC = (min)(A, c);
	ASSERT_TRUE( is_equal(AC, AC_r) );

	mat_t CB = (min)(c, B);
	ASSERT_TRUE( is_equal(CB, CB_r) );
}


MN_CASE( mat_arith_ex,  add_block_to_block )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	double v = 1.5;
	dblock<double> blk_x(LDim * n);
	dblock<double> blk_y(LDim * n, fill(v));

	for (index_t i = 0; i < LDim * n; ++i) blk_x[i] = double(i+1);

	cref_block<double, M, N> x(blk_x.ptr_data(), m, n, LDim);
	ref_block<double, M, N> y(blk_y.ptr_data(), m, n, LDim);

	dense_matrix<double, M, N> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r(i, j) = x(i, j) + v;
	}

	y += x;

	ASSERT_MAT_EQ(m, n, y, r);
}

MN_CASE( mat_arith_ex,  add_grid_to_grid )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	double v = 1.5;
	dblock<double> blk_x(LDim * n);
	dblock<double> blk_y(LDim * n, fill(v));

	for (index_t i = 0; i < LDim * n; ++i) blk_x[i] = double(i+1);

	cref_grid<double, M, N> x(blk_x.ptr_data(), m, n, rs, cs);
	ref_grid<double, M, N> y(blk_y.ptr_data(), m, n, rs, cs);

	dense_matrix<double, M, N> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r(i, j) = x(i, j) + v;
	}

	y += x;

	ASSERT_MAT_EQ(m, n, y, r);
}

MN_CASE( mat_arith_ex,  negate_block )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	double v = 1.5;
	dblock<double> blk_x(LDim * n);
	dblock<double> blk_y(LDim * n, fill(v));

	for (index_t i = 0; i < LDim * n; ++i) blk_x[i] = double(i+1);

	cref_block<double, M, N> x(blk_x.ptr_data(), m, n, LDim);
	ref_block<double, M, N> y(blk_y.ptr_data(), m, n, LDim);

	dense_matrix<double, M, N> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r(i, j) = - x(i, j);
	}

	y = -x;

	ASSERT_MAT_EQ(m, n, y, r);
}

MN_CASE( mat_arith_ex,  negate_grid )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	double v = 1.5;
	dblock<double> blk_x(LDim * n);
	dblock<double> blk_y(LDim * n, fill(v));

	for (index_t i = 0; i < LDim * n; ++i) blk_x[i] = double(i+1);

	cref_grid<double, M, N> x(blk_x.ptr_data(), m, n, rs, cs);
	ref_grid<double, M, N> y(blk_y.ptr_data(), m, n, rs, cs);

	dense_matrix<double, M, N> r(m, n);
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i) r(i, j) = - x(i, j);
	}

	y = -x;

	ASSERT_MAT_EQ(m, n, y, r);
}


BEGIN_TPACK( sch_choice_mat_mat )
	ADD_MN_CASE_3X3( sch_choice, mat_mat, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_mat_blk )
	ADD_MN_CASE_3X3( sch_choice, mat_blk, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_mat_grid )
	ADD_MN_CASE_3X3( sch_choice, mat_grid, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_blk_mat )
	ADD_MN_CASE_3X3( sch_choice, blk_mat, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_blk_blk )
	ADD_MN_CASE_3X3( sch_choice, blk_blk, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_blk_grid )
	ADD_MN_CASE_3X3( sch_choice, blk_grid, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_grid_mat )
	ADD_MN_CASE_3X3( sch_choice, grid_mat, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_grid_blk )
	ADD_MN_CASE_3X3( sch_choice, grid_blk, DM, DN )
END_TPACK

BEGIN_TPACK( sch_choice_grid_grid )
	ADD_MN_CASE_3X3( sch_choice, grid_grid, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_add )
	ADD_MN_CASE_3X3( mat_arith, add, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_add_ip )
	ADD_MN_CASE_3X3( mat_arith, add_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_sub )
	ADD_MN_CASE_3X3( mat_arith, sub, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_sub_ip )
	ADD_MN_CASE_3X3( mat_arith, sub_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_mul )
	ADD_MN_CASE_3X3( mat_arith, mul, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_mul_ip )
	ADD_MN_CASE_3X3( mat_arith, mul_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_div )
	ADD_MN_CASE_3X3( mat_arith, div, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_div_ip )
	ADD_MN_CASE_3X3( mat_arith, div_ip, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_neg )
	ADD_MN_CASE_3X3( mat_arith, neg, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_abs )
	ADD_MN_CASE_3X3( mat_arith, abs, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_max )
	ADD_MN_CASE_3X3( mat_arith, max, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_min )
	ADD_MN_CASE_3X3( mat_arith, min, DM, DN )
END_TPACK


BEGIN_TPACK( mat_arith_ex_add_block )
	ADD_MN_CASE_3X3( mat_arith_ex, add_block_to_block, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_ex_add_grid )
	ADD_MN_CASE_3X3( mat_arith_ex, add_grid_to_grid, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_ex_neg_block )
	ADD_MN_CASE_3X3( mat_arith_ex, negate_block, DM, DN )
END_TPACK

BEGIN_TPACK( mat_arith_ex_neg_grid )
	ADD_MN_CASE_3X3( mat_arith_ex, negate_grid, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( sch_choice_mat_mat )
	ADD_TPACK( sch_choice_mat_blk )
	ADD_TPACK( sch_choice_mat_grid )
	ADD_TPACK( sch_choice_blk_mat )
	ADD_TPACK( sch_choice_blk_blk )
	ADD_TPACK( sch_choice_blk_grid )
	ADD_TPACK( sch_choice_grid_mat )
	ADD_TPACK( sch_choice_grid_blk )
	ADD_TPACK( sch_choice_grid_grid )

	ADD_TPACK( mat_arith_add )
	ADD_TPACK( mat_arith_add_ip )
	ADD_TPACK( mat_arith_sub )
	ADD_TPACK( mat_arith_sub_ip )
	ADD_TPACK( mat_arith_mul )
	ADD_TPACK( mat_arith_mul_ip )
	ADD_TPACK( mat_arith_div )
	ADD_TPACK( mat_arith_div_ip )

	ADD_TPACK( mat_arith_neg )
	ADD_TPACK( mat_arith_abs )

	ADD_TPACK( mat_arith_max )
	ADD_TPACK( mat_arith_min )

	ADD_TPACK( mat_arith_ex_add_block )
	ADD_TPACK( mat_arith_ex_add_grid )
	ADD_TPACK( mat_arith_ex_neg_block )
	ADD_TPACK( mat_arith_ex_neg_grid )

END_MAIN_SUITE



