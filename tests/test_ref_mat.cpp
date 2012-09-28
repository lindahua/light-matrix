/**
 * @file test_ref_mat.cpp
 *
 * Unit testing of cref_matrix and ref_matrix
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;


// explicit instantiation

template class lmat::cref_matrix<double, 0, 0>;
template class lmat::cref_matrix<double, 0, 4>;
template class lmat::cref_matrix<double, 3, 0>;
template class lmat::cref_matrix<double, 3, 4>;

template class lmat::ref_matrix<double, 0, 0>;
template class lmat::ref_matrix<double, 0, 4>;
template class lmat::ref_matrix<double, 3, 0>;
template class lmat::ref_matrix<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT

static_assert(lmat::is_mat_xpr<lmat::cref_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::cref_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::cref_matrix<double> >::value, "Interface verification failed.");

static_assert(lmat::is_mat_xpr<lmat::ref_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::ref_matrix<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::ref_matrix<double> >::value, "Interface verification failed.");

#endif

MN_CASE( cref_mat, constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(m * n);
	const double *ps = s.ptr_data();

	cref_matrix<double, M, N> a(ps, m, n);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), m);
	ASSERT_EQ(a.size(), (size_t)a.nelems() );

	ASSERT_EQ(a.ptr_data(), ps);

	cref_matrix<double, M, N> a2(a);

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), m);
	ASSERT_EQ(a2.size(), (size_t)a2.nelems() );

	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( ref_mat, constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(m * n);
	double *ps = s.ptr_data();

	ref_matrix<double, M, N> a(ps, m, n);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), m);
	ASSERT_EQ(a.size(), (size_t)a.nelems() );

	ASSERT_EQ(a.ptr_data(), ps);

	ref_matrix<double, M, N> a2(a);

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), m);
	ASSERT_EQ(a2.size(), (size_t)a2.nelems() );

	ASSERT_EQ(a2.ptr_data(), ps);
}


MN_CASE( cref_mat, access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	cref_matrix<double, M, N> a(ref.ptr_data(), m, n);
	const cref_matrix<double, M, N>& ac = a;

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


MN_CASE( ref_mat, access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);

	for (index_t i = 0; i < m * n; ++i) ref[i] = double(i + 2);
	ref_matrix<double, M, N> a(ref.ptr_data(), m, n);
	const ref_matrix<double, M, N>& ac = a;

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


MN_CASE( ref_mat, assign )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s1(m * n);
	dblock<double> s2(m * n);

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < m * n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < m * n; ++i) s2[i] = double(2 * i + 3);

	ref_matrix<double, M, N> a1(ps1, m, n);
	ref_matrix<double, M, N> a2(ps2, m, n);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_VEC_EQ( m * n, a1, a2 );
}

MN_CASE( ref_mat, import )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(m * n);
	dblock<double> s(m * n, fill(-1.0));

	ref_matrix<double, M, N> a(s.ptr_data(), m, n);

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



BEGIN_TPACK( cref_mat_constructs )
	ADD_MN_CASE_3X3( cref_mat, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_constructs )
	ADD_MN_CASE_3X3( ref_mat, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( cref_mat_access )
	ADD_MN_CASE_3X3( cref_mat, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_access )
	ADD_MN_CASE_3X3( ref_mat, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_assign )
	ADD_MN_CASE_3X3( ref_mat, assign, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_import )
	ADD_MN_CASE_3X3( ref_mat, import, 3, 4 )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( cref_mat_constructs )
	ADD_TPACK( cref_mat_access )

	ADD_TPACK( ref_mat_constructs )
	ADD_TPACK( ref_mat_access )
	ADD_TPACK( ref_mat_assign )
	ADD_TPACK( ref_mat_import )
END_MAIN_SUITE




