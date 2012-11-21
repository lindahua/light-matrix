/**
 * @file test_ref_grid.cpp
 *
 * @brief Unit testing of classes ref_grid and cref_grid
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/ref_grid.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::cref_grid<double, 0, 0>;
template class lmat::cref_grid<double, 0, 4>;
template class lmat::cref_grid<double, 3, 0>;
template class lmat::cref_grid<double, 3, 4>;

template class lmat::ref_grid<double, 0, 0>;
template class lmat::ref_grid<double, 0, 4>;
template class lmat::ref_grid<double, 3, 0>;
template class lmat::ref_grid<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT

static_assert(lmat::is_mat_xpr<lmat::cref_grid<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::cref_grid<double> >::value, "Interface verification failed.");

static_assert(lmat::is_mat_xpr<lmat::ref_grid<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::ref_grid<double> >::value, "Interface verification failed.");

#endif


template<int M, int N>
inline void verify_layout(const cref_grid<double, M, N>& mat,
		index_t m, index_t n, index_t rs, index_t cs)
{
	ASSERT_EQ( mat.nrows(), m );
	ASSERT_EQ( mat.ncolumns(), n );
	ASSERT_EQ( mat.nelems(), m * n );
	ASSERT_EQ( mat.row_stride(), rs );
	ASSERT_EQ( mat.col_stride(), cs );
}

template<int M, int N>
inline void verify_layout(const ref_grid<double, M, N>& mat,
		index_t m, index_t n, index_t rs, index_t cs)
{
	ASSERT_EQ( mat.nrows(), m );
	ASSERT_EQ( mat.ncolumns(), n );
	ASSERT_EQ( mat.nelems(), m * n );
	ASSERT_EQ( mat.row_stride(), rs );
	ASSERT_EQ( mat.col_stride(), cs );
}


MN_CASE( cref_grid, constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;
	const index_t rs = 2;
	const index_t cs = 15;

	dblock<double> s(cs * n);
	const double *ps = s.ptr_data();

	cref_grid<double, M, N> a(ps, m, n, rs, cs);

	verify_layout(a, m, n, rs, cs);
	ASSERT_EQ(a.ptr_data(), ps);

	cref_grid<double, M, N> a2(a);
	verify_layout(a2, m, n, rs, cs);

	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( ref_grid, constructs )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;
	const index_t rs = 2;
	const index_t cs = 15;

	dblock<double> s(cs * n);
	double *ps = s.ptr_data();

	ref_grid<double, M, N> a(ps, m, n, rs, cs);

	verify_layout(a, m, n, rs, cs);
	ASSERT_EQ(a.ptr_data(), ps);

	ref_grid<double, M, N> a2(a);

	verify_layout(a2, m, n, rs, cs);
	ASSERT_EQ(a2.ptr_data(), ps);
}


MN_CASE( cref_grid, access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;
	const index_t rs = 2;
	const index_t cs = 15;

	dblock<double> ref(cs * n);

	for (index_t i = 0; i < cs * n; ++i) ref[i] = double(i + 2);
	cref_grid<double, M, N> a(ref.ptr_data(), m, n, rs, cs);
	const cref_grid<double, M, N>& ac = a;

	for (index_t j = 0; j < n; ++j)
	{
		ASSERT_EQ(a.ptr_col(j), a.ptr_data() + j * cs);
		ASSERT_EQ(ac.ptr_col(j), ac.ptr_data() + j * cs);

		for (index_t i = 0; i < m; ++i)
		{
			double v0 = ref[i * rs + j * cs];

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
			ASSERT_EQ( a[i], ref[i * rs] );
			ASSERT_EQ( ac[i], ref[i * rs] );
		}
	}
	else if (M == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i * cs] );
			ASSERT_EQ( ac[i], ref[i * cs] );
		}
	}
}

MN_CASE( ref_grid, access )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;
	const index_t rs = 2;
	const index_t cs = 15;

	dblock<double> ref(cs * n);

	for (index_t i = 0; i < cs * n; ++i) ref[i] = double(i + 2);
	ref_grid<double, M, N> a(ref.ptr_data(), m, n, rs, cs);
	const ref_grid<double, M, N>& ac = a;

	for (index_t j = 0; j < n; ++j)
	{
		ASSERT_EQ(a.ptr_col(j), a.ptr_data() + j * cs);
		ASSERT_EQ(ac.ptr_col(j), ac.ptr_data() + j * cs);

		for (index_t i = 0; i < m; ++i)
		{
			double v0 = ref[i * rs + j * cs];

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
			ASSERT_EQ( a[i], ref[i * rs] );
			ASSERT_EQ( ac[i], ref[i * rs] );
		}
	}
	else if (M == 1)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ( a[i], ref[i * cs] );
			ASSERT_EQ( ac[i], ref[i * cs] );
		}
	}
}


MN_CASE( ref_grid, assign )
{

	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	const index_t rs1 = 2;
	const index_t rs2 = 3;
	const index_t cs1 = 15;
	const index_t cs2 = 17;

	dblock<double> s1(cs1 * n);
	dblock<double> s2(cs2 * n);

	double *ps1 = s1.ptr_data();
	double *ps2 = s2.ptr_data();

	for (index_t i = 0; i < cs1 * n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < cs2 * n; ++i) s2[i] = double(2 * i + 3);

	dblock<double> s1r(s1);

	ref_grid<double, M, N> a1(ps1, m, n, rs1, cs1);
	ref_grid<double, M, N> a2(ps2, m, n, rs2, cs2);

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
			s1r[i * rs1 + j * cs1] = s2[i * rs2 + j * cs2];
		}
	}

	ASSERT_VEC_EQ(cs1 * n, s1, s1r);
}


MN_CASE( ref_grid, import )
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;
	const index_t rs = 2;
	const index_t cs = 15;

	dblock<double> ref(cs * n);
	dblock<double> s0(cs * n, fill(-1.0));
	dblock<double> s(s0);
	double *ps = s.ptr_data();

	ref_grid<double, M, N> a(ps, m, n, rs, cs);

	// fill_value

	const double v1 = 2.5;
	a << v1;

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i * rs + j * cs] = v1;
	}

	verify_layout(a, m, n, rs, cs);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(cs * n, ps, ref);

	// copy_value

	dblock<double> csrc(m * n);

	for (index_t i = 0; i < m * n; ++i) csrc[i] = double(i + 2);
	a << csrc.ptr_data();

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i * rs + j * cs] = csrc[i + j * m];
	}

	verify_layout(a, m, n, rs, cs);
	ASSERT_EQ(a.ptr_data(), ps);
	ASSERT_VEC_EQ(cs * n, ps, ref);

}



BEGIN_TPACK( cref_grid_constructs )
	ADD_MN_CASE_3X3( cref_grid, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_grid_constructs )
	ADD_MN_CASE_3X3( ref_grid, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( cref_grid_access )
	ADD_MN_CASE_3X3( cref_grid, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_grid_access )
	ADD_MN_CASE_3X3( ref_grid, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_grid_assign )
	ADD_MN_CASE_3X3( ref_grid, assign, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_grid_import )
	ADD_MN_CASE_3X3( ref_grid, import, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( cref_grid_constructs )
	ADD_TPACK( cref_grid_access )

	ADD_TPACK( ref_grid_constructs )
	ADD_TPACK( ref_grid_access )
	ADD_TPACK( ref_grid_assign )
	ADD_TPACK( ref_grid_import )
END_MAIN_SUITE






