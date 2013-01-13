/**
 * @file test_step_vec.cpp
 *
 * @brief Unit testing of step vector classes
 *
 * @author Dahua Lin
 */


#include "../test_base.h"

#include <light_mat/matrix/step_vecs.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::cstep_col<double, 0>;
template class lmat::cstep_col<double, 4>;
template class lmat::cstep_row<double, 0>;
template class lmat::cstep_row<double, 4>;

template class lmat::step_col<double, 0>;
template class lmat::step_col<double, 4>;
template class lmat::step_row<double, 0>;
template class lmat::step_row<double, 4>;

template<class Vec>
inline void verify_steprow_layout(const Vec& a, index_t n, index_t step)
{
	ASSERT_EQ(a.nrows(), 1);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), n);
	ASSERT_EQ(a.row_stride(), 1);
	ASSERT_EQ(a.col_stride(), step);
}

template<class Vec>
inline void verify_stepcol_layout(const Vec& a, index_t m, index_t step)
{
	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), 1);
	ASSERT_EQ(a.nelems(), m);
	ASSERT_EQ(a.row_stride(), step);
	ASSERT_EQ(a.col_stride(), 0);
}


N_CASE( cstep_row_constructs )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s(n * step);
	const double *ps = s.ptr_data();

	cstep_row<double, N> a(ps, n, step);

	verify_steprow_layout(a, n, step);
	ASSERT_EQ(a.ptr_data(), ps);

	cstep_row<double, N> a2(a);

	verify_steprow_layout(a2, n, step);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( step_row_constructs )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s(n * step);
	double *ps = s.ptr_data();

	step_row<double, N> a(ps, n, step);

	verify_steprow_layout(a, n, step);
	ASSERT_EQ(a.ptr_data(), ps);

	step_row<double, N> a2(a);

	verify_steprow_layout(a2, n, step);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( cstep_col_constructs )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s(n * step);
	const double *ps = s.ptr_data();

	cstep_col<double, N> a(ps, n, step);

	verify_stepcol_layout(a, n, step);
	ASSERT_EQ(a.ptr_data(), ps);

	cstep_col<double, N> a2(a);

	verify_stepcol_layout(a2, n, step);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( step_col_constructs )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s(n * step);
	double *ps = s.ptr_data();

	step_col<double, N> a(ps, n, step);

	verify_stepcol_layout(a, n, step);
	ASSERT_EQ(a.ptr_data(), ps);

	step_col<double, N> a2(a);

	verify_stepcol_layout(a2, n, step);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( step_row_assign )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s1(n * step, zero());
	dblock<double> s2(n * step, zero());

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < n * step; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < n * step; ++i) s2[i] = double(2 * i + 3);

	step_row<double, N> a1(ps1, n, step);
	step_row<double, N> a2(ps2, n, step);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_VEC_EQ( n, a1, a2 );
}


N_CASE( step_col_assign )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t step = 2;

	dblock<double> s1(n * step, zero());
	dblock<double> s2(n * step, zero());

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < n * step; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < n * step; ++i) s2[i] = double(2 * i + 3);

	step_col<double, N> a1(ps1, n, step);
	step_col<double, N> a2(ps2, n, step);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_VEC_EQ( n, a1, a2 );
}


AUTO_TPACK( cstep_row_constructs )
{
	ADD_N_CASE_3( cstep_row_constructs, 4 )
}

AUTO_TPACK( step_row_constructs )
{
	ADD_N_CASE_3( step_row_constructs, 4 )
}

AUTO_TPACK( cstep_col_constructs )
{
	ADD_N_CASE_3( cstep_col_constructs, 4 )
}

AUTO_TPACK( step_col_constructs )
{
	ADD_N_CASE_3( step_col_constructs, 4 )
}


AUTO_TPACK( step_row_assign )
{
	ADD_N_CASE_3( step_row_assign, 4 )
}

AUTO_TPACK( step_col_assign )
{
	ADD_N_CASE_3( step_col_assign, 4 )
}





