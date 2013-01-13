/**
 * @file test_dense_mat.cpp
 *
 * Unit testing of dense matrix
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/dense_matrix.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::dense_matrix<double, 0, 0>;
template class lmat::dense_matrix<double, 0, 4>;
template class lmat::dense_matrix<double, 3, 0>;
template class lmat::dense_matrix<double, 3, 4>;

static_assert(lmat::meta::is_mat_xpr<lmat::dense_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::meta::is_regular_mat<lmat::dense_matrix<double> >::value, "Interface verification failed.");


template<int M, int N>
inline void verify_layout(const dense_matrix<double, M, N>& a, index_t m, index_t n)
{
	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.row_stride(), 1);
	ASSERT_EQ(a.col_stride(), m);
}


MN_CASE( dense_mat_constructs )
{
	// default construction

	dense_matrix<double, M, N> a0;

	verify_layout(a0, M, N);

	if (M > 0 && N > 0)
	{
		ASSERT_TRUE(a0.ptr_data() != 0);
	}
	else
	{
		ASSERT_TRUE(a0.ptr_data() == 0);
	}

	// size given construction

	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dense_matrix<double, M, N> a1(m, n);

	verify_layout(a1, m, n);
	ASSERT_TRUE(a1.ptr_data() != 0);

}


SIMPLE_CASE( dense_mat_initializes )
{
	const index_t m = 2;
	const index_t n = 3;

	dense_matrix<double> a(m, n, rm_({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}));

	double a0_src[m * n] = {1.0, 4.0, 2.0, 5.0, 3.0, 6.0};
	dense_matrix<double> a0(m, n, copy_from(a0_src));

	ASSERT_EQ( a.nrows(), m );
	ASSERT_EQ( a.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, a, a0 );


	dense_matrix<double> b(m, n, cm_({1.0, 2.0, 3.0, 4.0, 5.0, 6.0}));

	double b0_src[m * n] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0};
	dense_matrix<double> b0(m, n, copy_from(b0_src));

	ASSERT_EQ( b.nrows(), m );
	ASSERT_EQ( b.ncolumns(), n );
	ASSERT_MAT_EQ( m, n, b, b0 );
}


MN_CASE( dense_mat_generates )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	// zeros

	dense_matrix<double, M, N> a0(m, n, zero());
	for (index_t i = 0; i < m * n; ++i) ref[i] = double(0);

	verify_layout(a0, m, n);
	ASSERT_VEC_EQ(m * n, a0, ref);

	// fill_value

	const double v1 = 2.5;
	dense_matrix<double, M, N> a1(m, n, fill(v1));
	for (index_t i = 0; i < m * n; ++i) ref[i] = v1;

	verify_layout(a1, m, n);
	ASSERT_VEC_EQ(m * n, a1, ref);

	// copy_value

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> a2(m, n, copy_from(ref.ptr_data()));

	verify_layout(a2, m, n);
	ASSERT_VEC_EQ(m * n, a2, ref);
}


MN_CASE( dense_mat_copy_constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> a(m, n, copy_from(ref.ptr_data()));

	dense_matrix<double, M, N> a2(a);

	verify_layout(a2, m, n);
	ASSERT_NE(a.ptr_data(), a2.ptr_data());

	ASSERT_VEC_EQ(m * n, a, a2);
}


MN_CASE( dense_mat_move_constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dense_matrix<double, M, N> a(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = double(i + 2);

	double *p = a.ptr_data();

	dense_matrix<double, M, N> b = std::move(a);

	verify_layout(b, m, n);

	if (M == 0 || N == 0)
	{
		ASSERT_EQ( a.ptr_data(), nullptr );
		ASSERT_EQ( b.ptr_data(), p );
	}
	else
	{
		ASSERT_EQ( a.ptr_data(), p );
		ASSERT_MAT_EQ( m, n, a, b );
	}
}


MN_CASE( dense_mat_access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> a(m, n, copy_from(ref.ptr_data()));
	const dense_matrix<double, M, N>& ac = a;

	ASSERT_VEC_EQ(m * n, a, ref);

	for (index_t j = 0; j < n; ++j)
	{
		ASSERT_EQ(a.ptr_col(j), a.ptr_data() + j * m);
		ASSERT_EQ(ac.ptr_col(j), ac.ptr_data() + j * m);

		for (index_t i = 0; i < m; ++i)
		{
			double v0 = ref[i + j * m];

			ASSERT_EQ( a.elem(i, j), v0 );
			ASSERT_EQ( a(i, j), v0 );
			ASSERT_EQ( ac.elem(i, j), v0 );
			ASSERT_EQ( ac(i, j), v0 );
		}
	}

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ( a[i], ref[i] );
		ASSERT_EQ( ac[i], ref[i] );
	}
}

MN_CASE( dense_mat_resize )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const index_t m2 = M == 0 ? 4 : M;
	const index_t n2 = N == 0 ? 3 : N;

	const index_t m3 = M == 0 ? 5 : M;
	const index_t n3 = N == 0 ? 6 : N;

	dense_matrix<double, M, N> a(m, n);

	verify_layout(a, m, n);
	const double *p1 = a.ptr_data();

	a.require_size(m2, n2);

	verify_layout(a, m2, n2);
	const double *p2 = a.ptr_data();

	if (m2 * n2 == m * n)
	{
		ASSERT_EQ( p2, p1 );
	}
	else
	{
		ASSERT_NE( p2, p1 );
	}

	a.require_size(m3, n3);

	verify_layout(a, m3, n3);
	const double *p3 = a.ptr_data();

	if (m2 * n2 == m3 * n3)
	{
		ASSERT_EQ( p3, p2 );
	}
	else
	{
		ASSERT_NE( p3, p2 );
	}
}


MN_CASE( dense_mat_assign )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> s(m, n, copy_from(ref.ptr_data()));

	dense_matrix<double, M, N> a;

	a = s;

	verify_layout(a, m, n);
	ASSERT_NE( a.ptr_data(), 0 );
	ASSERT_NE( a.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( m * n, a, s );

	dense_matrix<double, M, N> b(m, n, zero());

	const double *pb = b.ptr_data();

	ASSERT_NE( pb, 0 );
	ASSERT_NE( pb, s.ptr_data() );

	b = s;

	verify_layout(b, m, n);
	ASSERT_EQ( b.ptr_data(), pb );


	ASSERT_VEC_EQ( m * n, b, s );

	const index_t m2 = M == 0 ? 5 : M;
	const index_t n2 = N == 0 ? 6 : N;

	dense_matrix<double, M, N> c(m2, n2, zero());

	c = s;

	verify_layout(c, m, n);
	ASSERT_NE( c.ptr_data(), 0 );
	ASSERT_NE( c.ptr_data(), s.ptr_data() );

	ASSERT_VEC_EQ( m * n, c, s );
}



MN_CASE( dense_mat_move_assign )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dense_matrix<double, M, N> a(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = double(i + 2);

	double *p = a.ptr_data();

	dense_matrix<double, M, N> b;
	b = std::move(a);

	verify_layout(b, m, n);

	if (M == 0 || N == 0)
	{
		ASSERT_EQ( a.ptr_data(), nullptr );
		ASSERT_EQ( b.ptr_data(), p );
	}
	else
	{
		ASSERT_EQ( a.ptr_data(), p );
		ASSERT_MAT_EQ( m, n, a, b );
	}

	dense_matrix<double, M, N> c(m, n);
	c = std::move(b);

	if (M == 0 || N == 0)
	{
		ASSERT_EQ( b.ptr_data(), nullptr );
		ASSERT_EQ( c.ptr_data(), p );
	}
	else
	{
		ASSERT_MAT_EQ( m, n, b, c );
	}

}


MN_CASE( dense_mat_import )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);
	dense_matrix<double, M, N> a(m, n, fill(-1.0));

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < m * n; ++i) ref[i] = v1;

	verify_layout(a, m, n);
	ASSERT_VEC_EQ(m * n, a, ref);

	// copy_value

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	verify_layout(a, m, n);
	ASSERT_VEC_EQ(m * n, a, ref);
}

MN_CASE( dense_mat_swap )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const index_t m2 = M == 0 ? 6 : M;
	const index_t n2 = N == 0 ? 5 : N;

	dblock<double> s(m * n);
	for (index_t i = 0; i < m * n; ++i) s[i] = double(i + 2);

	dblock<double> s2(m2 * n2);
	for (index_t i = 0; i < m2 * n2; ++i) s2[i] = double(2 * i + 3);

	dense_matrix<double, M, N> a(m, n, copy_from(s.ptr_data()));
	dense_matrix<double, M, N> a2(m2, n2, copy_from(s2.ptr_data()));

	const double *p = a.ptr_data();
	const double *p2 = a2.ptr_data();

	swap(a, a2);

	verify_layout(a, m2, n2);
	verify_layout(a2, m, n);

	if (M == 0 || N == 0)
	{
		ASSERT_EQ( a.ptr_data(), p2 );
		ASSERT_EQ( a2.ptr_data(), p );
	}
	else
	{
		ASSERT_EQ( a.ptr_data(), p );
		ASSERT_EQ( a2.ptr_data(), p2 );
	}

	ASSERT_VEC_EQ( m2 * n2, a, s2 );
	ASSERT_VEC_EQ( m * n, a2, s );
}

LTEST_INIT_AUTOSUITE

AUTO_TPACK( dense_mat_constructs )
{
	ADD_MN_CASE_3X3( dense_mat_constructs, 3, 4 )
	ADD_SIMPLE_CASE( dense_mat_initializes )
}

AUTO_TPACK( dense_mat_generates )
{
	ADD_MN_CASE_3X3( dense_mat_generates, 3, 4 )
}

AUTO_TPACK( dense_mat_copycon )
{
	ADD_MN_CASE_3X3( dense_mat_copy_constructs, 3, 4 )
}

AUTO_TPACK( dense_mat_movecon )
{
	ADD_MN_CASE_3X3( dense_mat_move_constructs, 3, 4 )
}

AUTO_TPACK( dense_mat_access )
{
	ADD_MN_CASE_3X3( dense_mat_access, 3, 4 )
}

AUTO_TPACK( dense_mat_resize )
{
	ADD_MN_CASE_3X3( dense_mat_resize, 3, 4 )
}

AUTO_TPACK( dense_mat_assign )
{
	ADD_MN_CASE_3X3( dense_mat_assign, 3, 4 )
}

AUTO_TPACK( dense_mat_move_assign )
{
	ADD_MN_CASE_3X3( dense_mat_move_assign, 3, 4 )
}

AUTO_TPACK( dense_mat_import )
{
	ADD_MN_CASE_3X3( dense_mat_import, 3, 4 )
}

AUTO_TPACK( dense_mat_swap )
{
	ADD_MN_CASE_3X3( dense_mat_swap, 3, 4 )
}

