/**
 * @file test_ref_block.cpp
 *
 * Unit testing of cref_block and ref_block
 *
 * @author Dahua Lin
 */


#include "test_base.h"

#include <light_mat/matrix/ref_block.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;


// explicit instantiation

template class lmat::cref_block<double, 0, 0>;
template class lmat::cref_block<double, 0, 4>;
template class lmat::cref_block<double, 3, 0>;
template class lmat::cref_block<double, 3, 4>;

template class lmat::ref_block<double, 0, 0>;
template class lmat::ref_block<double, 0, 4>;
template class lmat::ref_block<double, 3, 0>;
template class lmat::ref_block<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT

static_assert(lmat::meta::is_mat_xpr<lmat::cref_block<double> >::value, "Interface verification failed.");
static_assert(lmat::meta::is_dense_mat<lmat::cref_block<double> >::value, "Interface verification failed.");

static_assert(lmat::meta::is_mat_xpr<lmat::ref_block<double> >::value, "Interface verification failed.");
static_assert(lmat::meta::is_dense_mat<lmat::ref_block<double> >::value, "Interface verification failed.");

#endif


template<int M, int N>
inline void verify_layout(const cref_block<double, M, N>& mat,
		index_t m, index_t n, index_t ldim)
{
	ASSERT_EQ( mat.nrows(), m );
	ASSERT_EQ( mat.ncolumns(), n );
	ASSERT_EQ( mat.nelems(), m * n );
	ASSERT_EQ( mat.row_stride(), 1 );
	ASSERT_EQ( mat.col_stride(), ldim );
}

template<int M, int N>
inline void verify_layout(const ref_block<double, M, N>& mat,
		index_t m, index_t n, index_t ldim)
{
	ASSERT_EQ( mat.nrows(), m );
	ASSERT_EQ( mat.ncolumns(), n );
	ASSERT_EQ( mat.nelems(), m * n );
	ASSERT_EQ( mat.row_stride(), 1 );
	ASSERT_EQ( mat.col_stride(), ldim );
}


MN_CASE( cref_block, constructs )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(ldim * n);
	const double *ps = s.ptr_data();

	cref_block<double, M, N> a(ps, m, n, ldim);

	verify_layout(a, m, n, ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	cref_block<double, M, N> a2(a);
	verify_layout(a2, m, n, ldim);

	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( ref_block, constructs )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s(ldim * n);
	double *ps = s.ptr_data();

	ref_block<double, M, N> a(ps, m, n, ldim);

	verify_layout(a, m, n, ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ref_block<double, M, N> a2(a);

	verify_layout(a2, m, n, ldim);
	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( cref_block, access )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(ldim * n);

	for (index_t i = 0; i < ldim * n; ++i) ref[i] = double(i + 2);
	cref_block<double, M, N> a(ref.ptr_data(), m, n, ldim);
	const cref_block<double, M, N>& ac = a;

	for (index_t j = 0; j < n; ++j)
	{
		ASSERT_EQ(a.ptr_col(j), a.ptr_data() + j * ldim);
		ASSERT_EQ(ac.ptr_col(j), ac.ptr_data() + j * ldim);

		for (index_t i = 0; i < m; ++i)
		{
			double v0 = ref[i + j * ldim];

			ASSERT_EQ( a.elem(i, j), v0 );
			ASSERT_EQ( a(i, j), v0 );
			ASSERT_EQ( ac.elem(i, j), v0 );
			ASSERT_EQ( ac(i, j), v0 );
		}
	}

	if (N == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i] );
			ASSERT_EQ( ac[i], ref[i] );
		}
	}
	else if (M == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i * ldim] );
			ASSERT_EQ( ac[i], ref[i * ldim] );
		}
	}
}


MN_CASE( ref_block, access )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(ldim * n);

	for (index_t i = 0; i < ldim * n; ++i) ref[i] = double(i + 2);
	ref_block<double, M, N> a(ref.ptr_data(), m, n, ldim);
	const ref_block<double, M, N>& ac = a;

	for (index_t j = 0; j < n; ++j)
	{
		ASSERT_EQ(a.ptr_col(j), a.ptr_data() + j * ldim);
		ASSERT_EQ(ac.ptr_col(j), ac.ptr_data() + j * ldim);

		for (index_t i = 0; i < m; ++i)
		{
			double v0 = ref[i + j * ldim];

			ASSERT_EQ( a.elem(i, j), v0 );
			ASSERT_EQ( a(i, j), v0 );
			ASSERT_EQ( ac.elem(i, j), v0 );
			ASSERT_EQ( ac(i, j), v0 );
		}
	}

	if (N == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i] );
			ASSERT_EQ( ac[i], ref[i] );
		}
	}
	else if (M == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i * ldim] );
			ASSERT_EQ( ac[i], ref[i * ldim] );
		}
	}
}


MN_CASE( ref_block, assign )
{
	const index_t ldim1 = 7;
	const index_t ldim2 = 9;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> s1(ldim1 * n);
	dblock<double> s2(ldim2 * n);

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < ldim1 * n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < ldim2 * n; ++i) s2[i] = double(2 * i + 3);

	dblock<double> s1r(s1);

	ref_block<double, M, N> a1(ps1, m, n, ldim1);
	ref_block<double, M, N> a2(ps2, m, n, ldim2);

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );
	ASSERT_NE( ps1, ps2 );

	a1 = a2;

	ASSERT_EQ( a1.ptr_data(), ps1 );
	ASSERT_EQ( a2.ptr_data(), ps2 );

	ASSERT_MAT_EQ(m, n, a1, a2);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			s1r[i + j * ldim1] = s2[i + j * ldim2];
		}
	}

	ASSERT_VEC_EQ(ldim1 * n, s1, s1r);
}


MN_CASE( ref_block, import )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> ref(ldim * n);
	dblock<double> s0(ldim * n, fill(-1.0));
	dblock<double> s(s0);
	double *ps = s.ptr_data();

	ref_block<double, M, N> a(ps, m, n, ldim);

	// fill_value

	const double v1 = 2.5;
	a << v1;

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim] = v1;
	}

	verify_layout(a, m, n, ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(ldim * n, ps, ref);

	// copy_value

	dblock<double> cs(m * n);

	for (index_t i = 0; i < m * n; ++i) cs[i] = double(i + 2);
	a << cs.ptr_data();

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim] = cs[i + j * m];
	}

	verify_layout(a, m, n, ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(ldim * n, ps, ref);

}



BEGIN_TPACK( cref_block_constructs )
	ADD_MN_CASE_3X3( cref_block, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_block_constructs )
	ADD_MN_CASE_3X3( ref_block, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( cref_block_access )
	ADD_MN_CASE_3X3( cref_block, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_block_access )
	ADD_MN_CASE_3X3( ref_block, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_block_assign )
	ADD_MN_CASE_3X3( ref_block, assign, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_block_import )
	ADD_MN_CASE_3X3( ref_block, import, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( cref_block_constructs )
	ADD_TPACK( cref_block_access )

	ADD_TPACK( ref_block_constructs )
	ADD_TPACK( ref_block_access )
	ADD_TPACK( ref_block_assign )
	ADD_TPACK( ref_block_import )
END_MAIN_SUITE


