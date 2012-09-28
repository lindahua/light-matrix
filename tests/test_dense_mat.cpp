/**
 * @file test_dense_mat.cpp
 *
 * Unit testing of dense matrix
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/dense_matrix.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::dense_matrix<double, 0, 0>;
template class lmat::dense_matrix<double, 0, 4>;
template class lmat::dense_matrix<double, 3, 0>;
template class lmat::dense_matrix<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT

static_assert(lmat::is_mat_xpr<lmat::dense_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::dense_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::dense_matrix<double> >::value, "Interface verification failed.");

#endif


MN_CASE( dense_mat, constructs )
{
	// default construction

	dense_matrix<double, M, N> a0;

	ASSERT_EQ(a0.nrows(), M);
	ASSERT_EQ(a0.ncolumns(), N);
	ASSERT_EQ(a0.nelems(), M * N);
	ASSERT_EQ(a0.lead_dim(), M);

	ASSERT_EQ(a0.size(), (size_t)a0.nelems() );

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

	ASSERT_EQ(a1.nrows(), m);
	ASSERT_EQ(a1.ncolumns(), n);
	ASSERT_EQ(a1.nelems(), m * n);
	ASSERT_EQ(a1.lead_dim(), m);

	ASSERT_EQ(a1.size(), (size_t)a1.nelems() );
	ASSERT_TRUE(a1.ptr_data() != 0);

}


MN_CASE( dense_mat, generates )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	// zeros

	dense_matrix<double, M, N> a0(m, n, zero());
	for (index_t i = 0; i < m * n; ++i) ref[i] = double(0);

	ASSERT_EQ(a0.nrows(), m);
	ASSERT_EQ(a0.ncolumns(), n);
	ASSERT_EQ(a0.nelems(), m * n);
	ASSERT_EQ(a0.lead_dim(), m);
	ASSERT_VEC_EQ(m * n, a0, ref);

	// fill_value

	const double v1 = 2.5;
	dense_matrix<double, M, N> a1(m, n, fill(v1));
	for (index_t i = 0; i < m * n; ++i) ref[i] = v1;

	ASSERT_EQ(a1.nrows(), m);
	ASSERT_EQ(a1.ncolumns(), n);
	ASSERT_EQ(a1.nelems(), m * n);
	ASSERT_EQ(a1.lead_dim(), m);
	ASSERT_VEC_EQ(m * n, a1, ref);

	// copy_value

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> a2(m, n, copy_from(ref.ptr_data()));

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), m);
	ASSERT_VEC_EQ(m * n, a2, ref);
}


MN_CASE( dense_mat, copy_constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> a(m, n, copy_from(ref.ptr_data()));

	dense_matrix<double, M, N> a2(a);

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), m);

	ASSERT_NE(a.ptr_data(), a2.ptr_data());

	ASSERT_VEC_EQ(m * n, a, a2);
}


MN_CASE( dense_mat, access )
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

MN_CASE( dense_mat, resize )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const index_t m2 = M == 0 ? 4 : M;
	const index_t n2 = N == 0 ? 3 : N;

	const index_t m3 = M == 0 ? 5 : M;
	const index_t n3 = N == 0 ? 6 : N;

	dense_matrix<double, M, N> a(m, n);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), m);

	const double *p1 = a.ptr_data();

	a.require_size(m2, n2);

	ASSERT_EQ(a.nrows(), m2);
	ASSERT_EQ(a.ncolumns(), n2);
	ASSERT_EQ(a.nelems(), m2 * n2);
	ASSERT_EQ(a.lead_dim(), m2);

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

	ASSERT_EQ(a.nrows(), m3);
	ASSERT_EQ(a.ncolumns(), n3);
	ASSERT_EQ(a.nelems(), m3 * n3);
	ASSERT_EQ(a.lead_dim(), m3);

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


MN_CASE( dense_mat, assign )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	dense_matrix<double, M, N> s(m, n, copy_from(ref.ptr_data()));

	dense_matrix<double, M, N> a;

	a = s;

	ASSERT_NE( a.ptr_data(), 0 );
	ASSERT_NE( a.ptr_data(), s.ptr_data() );
	ASSERT_EQ( a.nrows(), m);
	ASSERT_EQ( a.ncolumns(), n);
	ASSERT_EQ( a.nelems(), m * n);
	ASSERT_EQ( a.lead_dim(), m);

	ASSERT_VEC_EQ( m * n, a, s );

	dense_matrix<double, M, N> b(m, n, zero());

	const double *pb = b.ptr_data();

	ASSERT_NE( pb, 0 );
	ASSERT_NE( pb, s.ptr_data() );

	b = s;

	ASSERT_EQ( b.ptr_data(), pb );
	ASSERT_EQ( b.nrows(), m);
	ASSERT_EQ( b.ncolumns(), n);
	ASSERT_EQ( b.nelems(), m * n);
	ASSERT_EQ( b.lead_dim(), m);

	ASSERT_VEC_EQ( m * n, b, s );

	const index_t m2 = M == 0 ? 5 : M;
	const index_t n2 = N == 0 ? 6 : N;

	dense_matrix<double, M, N> c(m2, n2, zero());

	c = s;

	ASSERT_NE( c.ptr_data(), 0 );
	ASSERT_NE( c.ptr_data(), s.ptr_data() );
	ASSERT_EQ( c.nrows(), m);
	ASSERT_EQ( c.ncolumns(), n);
	ASSERT_EQ( c.nelems(), m * n);
	ASSERT_EQ( c.lead_dim(), m);

	ASSERT_VEC_EQ( m * n, c, s );
}


MN_CASE( dense_mat, import )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);
	dense_matrix<double, M, N> a(m, n, fill(-1.0));

	// fill_value

	const double v1 = 2.5;
	a << v1;
	for (index_t i = 0; i < m * n; ++i) ref[i] = v1;

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), m);
	ASSERT_VEC_EQ(m * n, a, ref);

	// copy_value

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	a << ref.ptr_data();

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), m);
	ASSERT_VEC_EQ(m * n, a, ref);
}

MN_CASE( dense_mat, swap )
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

	ASSERT_EQ( a.nrows(), m2 );
	ASSERT_EQ( a.ncolumns(), n2 );
	ASSERT_EQ( a.nelems(), m2 * n2 );

	ASSERT_EQ( a2.nrows(), m );
	ASSERT_EQ( a2.ncolumns(), n );
	ASSERT_EQ( a2.nelems(), m * n );

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


BEGIN_TPACK( dense_mat_constructs )
	ADD_MN_CASE_3X3( dense_mat, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_generates )
	ADD_MN_CASE_3X3( dense_mat, generates, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_copycon )
	ADD_MN_CASE_3X3( dense_mat, copy_constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_access )
	ADD_MN_CASE_3X3( dense_mat, access, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_resize )
	ADD_MN_CASE_3X3( dense_mat, resize, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_assign )
	ADD_MN_CASE_3X3( dense_mat, assign, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_import )
	ADD_MN_CASE_3X3( dense_mat, import, 3, 4 )
END_TPACK

BEGIN_TPACK( dense_mat_swap )
	ADD_MN_CASE_3X3( dense_mat, swap, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( dense_mat_constructs )
	ADD_TPACK( dense_mat_generates )
	ADD_TPACK( dense_mat_copycon )
	ADD_TPACK( dense_mat_access )
	ADD_TPACK( dense_mat_resize )
	ADD_TPACK( dense_mat_assign )
	ADD_TPACK( dense_mat_import )
	ADD_TPACK( dense_mat_swap )
END_MAIN_SUITE


