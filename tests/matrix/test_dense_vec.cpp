/**
 * @file test_dense_vec.cpp
 *
 * Unit testing of dense_col and dense_row
 *
 * @author Dahua Lin
 */


#include "../test_base.h"

#include <light_mat/matrix/dense_matrix.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::dense_col<double, 0>;
template class lmat::dense_col<double, 4>;
template class lmat::dense_row<double, 0>;
template class lmat::dense_row<double, 4>;

static_assert(std::is_base_of<
		lmat::dense_matrix<double, 0, 1>,
		lmat::dense_col<double, 0> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::dense_matrix<double, 4, 1>,
		lmat::dense_col<double, 4> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::dense_matrix<double, 1, 0>,
		lmat::dense_row<double, 0> >::value, "Base verification failed.");

static_assert(std::is_base_of<
		lmat::dense_matrix<double, 1, 4>,
		lmat::dense_row<double, 4> >::value, "Base verification failed.");


template<int M, int N>
inline void verify_layout(const dense_matrix<double, M, N>& a, index_t m, index_t n)
{
	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.row_stride(), 1);
	ASSERT_EQ(a.col_stride(), m);
}


N_CASE( dense_col_constructs )
{
	// default construction

	dense_col<double, N> a0;

	verify_layout(a0, N, 1);
	if (N > 0)
	{
		ASSERT_TRUE(a0.ptr_data() != 0);
	}
	else
	{
		ASSERT_TRUE(a0.ptr_data() == 0);
	}

	// size given construction

	const index_t n = N == 0 ? 4 : N;

	dense_col<double, N> a1(n);
	verify_layout(a1, n, 1);

	ASSERT_TRUE(a1.ptr_data() != 0);
}


SIMPLE_CASE( dense_col_initializes )
{
	dense_col<double> a {1.0, 2.0, 3.0, 4.0, 5.0};
	double s[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

	ASSERT_EQ( a.nrows(), 5 );
	ASSERT_EQ( a.ncolumns(), 1 );
	ASSERT_VEC_EQ( 5, a, s );
}


N_CASE( dense_col_generates )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	// zeros

	dense_col<double, N> a0(n, zero());
	for (index_t i = 0; i < n; ++i) ref[i] = double(0);

	verify_layout(a0, n, 1);
	ASSERT_VEC_EQ(n, a0, ref);

	// fill_value

	const double v1 = 2.5;
	dense_col<double, N> a1(n, fill(v1));
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	verify_layout(a1, n, 1);
	ASSERT_VEC_EQ(n, a1, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_col<double, N> a2(n, copy_from(ref.ptr_data()));

	verify_layout(a2, n, 1);
	ASSERT_VEC_EQ(n, a2, ref);
}

N_CASE( dense_col_copy_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_col<double, N> a(n, copy_from(ref.ptr_data()));

	dense_col<double, N> a2(a);

	verify_layout(a2, n, 1);
	ASSERT_NE(a.ptr_data(), a2.ptr_data());

	ASSERT_VEC_EQ(n, a, a2);
}


N_CASE( dense_col_resize )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t n2 = N == 0 ? 5 : N;

	dense_col<double, N> a(n);
	const double *p1 = a.ptr_data();

	a.require_size(n);

	verify_layout(a, n, 1);
	ASSERT_EQ(a.ptr_data(), p1);

	a.require_size(n2);

	verify_layout(a, n2, 1);
	const double *p2 = a.ptr_data();

	if (n2 == n)
	{
		ASSERT_EQ( p2, p1 );
	}
	else
	{
		ASSERT_NE( p2, p1 );
	}
}


N_CASE( dense_col_assign )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_col<double, N> s(n, copy_from(ref.ptr_data()));

	dense_col<double, N> a;

	a = s;

	verify_layout(a, n, 1);
	ASSERT_NE( a.ptr_data(), 0 );
	ASSERT_NE( a.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( n, a, s );

	dense_col<double, N> b(n, zero());

	const double *pb = b.ptr_data();

	ASSERT_NE( pb, 0 );
	ASSERT_NE( pb, s.ptr_data() );

	b = s;

	verify_layout(b, n, 1);
	ASSERT_EQ( b.ptr_data(), pb );

	ASSERT_VEC_EQ( n, b, s );

	const index_t n2 = N == 0 ? 6 : N;

	dense_col<double, N> c(n2, zero());

	c = s;

	verify_layout(c, n, 1);
	ASSERT_NE( c.ptr_data(), 0 );
	ASSERT_NE( c.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( n, c, s );
}


N_CASE( dense_col_import )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	dense_col<double, N> a(n, fill(-1.0));

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	verify_layout(a, n, 1);
	ASSERT_VEC_EQ(n, a, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	verify_layout(a, n, 1);
	ASSERT_VEC_EQ(n, a, ref);
}

N_CASE( dense_col_swap )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t n2 = N == 0 ? 5 : N;

	dblock<double> s(n);
	for (index_t i = 0; i < n; ++i) s[i] = double(i + 2);

	dblock<double> s2(n2);
	for (index_t i = 0; i < n2; ++i) s2[i] = double(2 * i + 3);

	dense_col<double, N> a(n, copy_from(s.ptr_data()));
	dense_col<double, N> a2(n2, copy_from(s2.ptr_data()));

	const double *p = a.ptr_data();
	const double *p2 = a2.ptr_data();

	swap(a, a2);

	verify_layout(a, n2, 1);
	verify_layout(a2, n, 1);

	if (N == 0)
	{
		ASSERT_EQ( a.ptr_data(), p2 );
		ASSERT_EQ( a2.ptr_data(), p );
	}
	else
	{
		ASSERT_EQ( a.ptr_data(), p );
		ASSERT_EQ( a2.ptr_data(), p2 );
	}

	ASSERT_VEC_EQ( n2, a, s2 );
	ASSERT_VEC_EQ( n, a2, s );
}



N_CASE( dense_row_constructs )
{
	// default construction

	dense_row<double, N> a0;

	verify_layout(a0, 1, N);
	if (N > 0)
	{
		ASSERT_TRUE(a0.ptr_data() != 0);
	}
	else
	{
		ASSERT_TRUE(a0.ptr_data() == 0);
	}

	// size given construction

	const index_t n = N == 0 ? 4 : N;

	dense_row<double, N> a1(n);

	verify_layout(a1, 1, n);
	ASSERT_TRUE(a1.ptr_data() != 0);
}

SIMPLE_CASE( dense_row_initializes )
{
	dense_row<double> a {1.0, 2.0, 3.0, 4.0, 5.0};
	double s[5] = {1.0, 2.0, 3.0, 4.0, 5.0};

	ASSERT_EQ( a.nrows(), 1 );
	ASSERT_EQ( a.ncolumns(), 5 );
	ASSERT_VEC_EQ( 5, a, s );
}


N_CASE( dense_row_generates )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	// zeros

	dense_row<double, N> a0(n, zero());
	for (index_t i = 0; i < n; ++i) ref[i] = double(0);

	verify_layout(a0, 1, n);
	ASSERT_VEC_EQ(n, a0, ref);

	// fill_value

	const double v1 = 2.5;
	dense_row<double, N> a1(n, fill(v1));
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	verify_layout(a1, 1, n);
	ASSERT_VEC_EQ(n, a1, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_row<double, N> a2(n, copy_from(ref.ptr_data()));

	verify_layout(a2, 1, n);
	ASSERT_VEC_EQ(n, a2, ref);
}

N_CASE( dense_row_copy_constructs )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_row<double, N> a(n, copy_from(ref.ptr_data()));

	dense_row<double, N> a2(a);

	verify_layout(a2, 1, n);
	ASSERT_NE(a.ptr_data(), a2.ptr_data());

	ASSERT_VEC_EQ(n, a, a2);
}

N_CASE( dense_row_resize )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t n2 = N == 0 ? 5 : N;

	dense_row<double, N> a(n);
	const double *p1 = a.ptr_data();

	a.require_size(n);

	verify_layout(a, 1, n);
	ASSERT_EQ(a.ptr_data(), p1);

	a.require_size(n2);

	verify_layout(a, 1, n2);
	const double *p2 = a.ptr_data();

	if (n2 == n)
	{
		ASSERT_EQ( p2, p1 );
	}
	else
	{
		ASSERT_NE( p2, p1 );
	}

}


N_CASE( dense_row_assign )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	dense_row<double, N> s(n, copy_from(ref.ptr_data()));

	dense_row<double, N> a;

	a = s;

	verify_layout(a, 1, n);
	ASSERT_NE( a.ptr_data(), 0 );
	ASSERT_NE( a.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( n, a, s );

	dense_row<double, N> b(n, zero());

	const double *pb = b.ptr_data();

	ASSERT_NE( pb, 0 );
	ASSERT_NE( pb, s.ptr_data() );

	b = s;

	verify_layout(b, 1, n);
	ASSERT_EQ( b.ptr_data(), pb );

	ASSERT_VEC_EQ( n, b, s );

	const index_t n2 = N == 0 ? 6 : N;

	dense_row<double, N> c(n2, zero());

	c = s;

	verify_layout(c, 1, n);
	ASSERT_NE( c.ptr_data(), 0 );
	ASSERT_NE( c.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( n, c, s );
}

N_CASE( dense_row_import )
{
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(n);
	dense_row<double, N> a(n, fill(-1.0));

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < n; ++i) ref[i] = v1;

	verify_layout(a, 1, n);
	ASSERT_VEC_EQ(n, a, ref);

	// copy_value

	for (index_t i = 0; i < n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	verify_layout(a, 1, n);
	ASSERT_VEC_EQ(n, a, ref);
}

N_CASE( dense_row_swap )
{
	const index_t n = N == 0 ? 4 : N;
	const index_t n2 = N == 0 ? 5 : N;

	dblock<double> s(n);
	for (index_t i = 0; i < n; ++i) s[i] = double(i + 2);

	dblock<double> s2(n2);
	for (index_t i = 0; i < n2; ++i) s2[i] = double(2 * i + 3);

	dense_row<double, N> a(n, copy_from(s.ptr_data()));
	dense_row<double, N> a2(n2, copy_from(s2.ptr_data()));

	const double *p = a.ptr_data();
	const double *p2 = a2.ptr_data();

	swap(a, a2);

	verify_layout(a, 1, n2);
	verify_layout(a2, 1, n);

	if (N == 0)
	{
		ASSERT_EQ( a.ptr_data(), p2 );
		ASSERT_EQ( a2.ptr_data(), p );
	}
	else
	{
		ASSERT_EQ( a.ptr_data(), p );
		ASSERT_EQ( a2.ptr_data(), p2 );
	}

	ASSERT_VEC_EQ( n2, a, s2 );
	ASSERT_VEC_EQ( n, a2, s );
}




AUTO_TPACK( dense_col_constructs )
{
	ADD_N_CASE_3( dense_col_constructs, 4 )
	ADD_SIMPLE_CASE( dense_col_initializes )
}

AUTO_TPACK( dense_col_generates )
{
	ADD_N_CASE_3( dense_col_generates, 4 )
}

AUTO_TPACK( dense_col_copycon )
{
	ADD_N_CASE_3( dense_col_copy_constructs, 4 )
}

AUTO_TPACK( dense_col_resize )
{
	ADD_N_CASE_3( dense_col_resize, 4 )
}

AUTO_TPACK( dense_col_assign )
{
	ADD_N_CASE_3( dense_col_assign, 4 )
}

AUTO_TPACK( dense_col_import )
{
	ADD_N_CASE_3( dense_col_import, 4 )
}

AUTO_TPACK( dense_col_swap )
{
	ADD_N_CASE_3( dense_col_swap, 4 )
}


AUTO_TPACK( dense_row_constructs )
{
	ADD_N_CASE_3( dense_row_constructs, 4 )
}

AUTO_TPACK( dense_row_generates )
{
	ADD_N_CASE_3( dense_row_generates, 4 )
}

AUTO_TPACK( dense_row_copycon )
{
	ADD_N_CASE_3( dense_row_copy_constructs, 4 )
}

AUTO_TPACK( dense_row_resize )
{
	ADD_N_CASE_3( dense_row_resize, 4 )
}

AUTO_TPACK( dense_row_assign )
{
	ADD_N_CASE_3( dense_row_assign, 4 )
}

AUTO_TPACK( dense_row_import )
{
	ADD_N_CASE_3( dense_row_import, 4 )
}

AUTO_TPACK( dense_row_swap )
{
	ADD_N_CASE_3( dense_row_swap, 4 )
}




