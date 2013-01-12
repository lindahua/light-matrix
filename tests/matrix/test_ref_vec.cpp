/**
 * @file test_ref_col.cpp
 *
 * Unit testing of cref_col and ref_col class.
 *
 * @author Dahua Lin
 */


#include "../test_base.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::cref_col<double, 0>;
template class lmat::cref_col<double, 4>;
template class lmat::cref_row<double, 0>;
template class lmat::cref_row<double, 4>;

template class lmat::ref_col<double, 0>;
template class lmat::ref_col<double, 4>;
template class lmat::ref_row<double, 0>;
template class lmat::ref_row<double, 4>;


static_assert(std::is_base_of<
		lmat::ref_matrix<double, 0, 1>,
		lmat::ref_col<double, 0> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::ref_matrix<double, 4, 1>,
		lmat::ref_col<double, 4> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::ref_matrix<double, 1, 0>,
		lmat::ref_row<double, 0> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::ref_matrix<double, 1, 4>,
		lmat::ref_row<double, 4> >::value, "Base verification failed.");


template<int M, int N>
inline void verify_layout(const cref_matrix<double, M, N>& a, index_t m, index_t n)
{
	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.row_stride(), 1);
	ASSERT_EQ(a.col_stride(), m);
}

template<int M, int N>
inline void verify_layout(const ref_matrix<double, M, N>& a, index_t m, index_t n)
{
	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.row_stride(), 1);
	ASSERT_EQ(a.col_stride(), m);
}


N_CASE( cref_col_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(n);
	const double *ps = s.ptr_data();

	cref_col<double, N> a(ps, n);

	verify_layout(a, n, 1);
	ASSERT_EQ(a.ptr_data(), ps);

	cref_col<double, N> a2(a);

	verify_layout(a2, n, 1);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( ref_col_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(n);
	double *ps = s.ptr_data();

	ref_col<double, N> a(ps, n);

	verify_layout(a, n, 1);
	ASSERT_EQ(a.ptr_data(), ps);

	ref_col<double, N> a2(a);

	verify_layout(a2, n, 1);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( ref_col_assign )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s1(n);
	dblock<double> s2(n);

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < n; ++i) s2[i] = double(2 * i + 3);

	ref_col<double, N> a1(ps1, n);
	ref_col<double, N> a2(ps2, n);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_VEC_EQ( n, a1, a2 );
}


N_CASE( ref_col_import )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);
	dblock<double> s(n, fill(-1.0));

	double *ps = s.ptr_data();
	ref_col<double, N> a(ps, n);

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	ASSERT_EQ(a.ptr_data(), ps);
	ASSERT_VEC_EQ(n, a, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	ASSERT_EQ(a.ptr_data(), ps);
	ASSERT_VEC_EQ(n, a, ref);
}


N_CASE( cref_row_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(n);
	const double *ps = s.ptr_data();

	cref_row<double, N> a(ps, n);

	verify_layout(a, 1, n);
	ASSERT_EQ(a.ptr_data(), ps);

	cref_row<double, N> a2(a);

	verify_layout(a2, 1, n);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( ref_row_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(n);
	double *ps = s.ptr_data();

	ref_row<double, N> a(ps, n);

	verify_layout(a, 1, n);
	ASSERT_EQ(a.ptr_data(), ps);

	ref_row<double, N> a2(a);

	verify_layout(a2, 1, n);
	ASSERT_EQ(a.ptr_data(), ps);
}


N_CASE( ref_row_assign )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s1(n);
	dblock<double> s2(n);

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < n; ++i) s2[i] = double(2 * i + 3);

	ref_row<double, N> a1(ps1, n);
	ref_row<double, N> a2(ps2, n);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_VEC_EQ( n, a1, a2 );
}


N_CASE( ref_row_import )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);
	dblock<double> s(n, fill(-1.0));

	double *ps = s.ptr_data();
	ref_row<double, N> a(ps, n);

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	ASSERT_EQ(a.ptr_data(), ps);
	ASSERT_VEC_EQ(n, a, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	ASSERT_EQ(a.ptr_data(), ps);
	ASSERT_VEC_EQ(n, a, ref);
}



AUTO_TPACK( cref_col_constructs )
{
	ADD_N_CASE_3( cref_col_constructs, 4 )
}

AUTO_TPACK( ref_col_constructs )
{
	ADD_N_CASE_3( ref_col_constructs, 4 )
}

AUTO_TPACK( ref_col_assign )
{
	ADD_N_CASE_3( ref_col_assign, 4 )
}

AUTO_TPACK( ref_col_import )
{
	ADD_N_CASE_3( ref_col_import, 4 )
}


AUTO_TPACK( cref_row_constructs )
{
	ADD_N_CASE_3( cref_row_constructs, 4 )
}

AUTO_TPACK( ref_row_constructs )
{
	ADD_N_CASE_3( ref_row_constructs, 4 )
}

AUTO_TPACK( ref_row_assign )
{
	ADD_N_CASE_3( ref_row_assign, 4 )
}

AUTO_TPACK( ref_row_import )
{
	ADD_N_CASE_3( ref_row_import, 4 )
}



