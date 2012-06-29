/**
 * @file test_ref_mat_ex.cpp
 *
 * Unit testing of cref_matrix_ex and ref_matrix_ex
 *
 * @author Dahua Lin
 */


#include "test_base.h"

#include <light_mat/matrix/ref_matrix_ex.h>
#include <light_mat/core/array.h>

using namespace lmat;
using namespace lmat::test;


// explicit instantiation

template class lmat::cref_matrix_ex<double, DynamicDim, DynamicDim>;
template class lmat::cref_matrix_ex<double, DynamicDim, 4>;
template class lmat::cref_matrix_ex<double, 3, DynamicDim>;
template class lmat::cref_matrix_ex<double, 3, 4>;

template class lmat::ref_matrix_ex<double, DynamicDim, DynamicDim>;
template class lmat::ref_matrix_ex<double, DynamicDim, 4>;
template class lmat::ref_matrix_ex<double, 3, DynamicDim>;
template class lmat::ref_matrix_ex<double, 3, 4>;

#ifdef LMAT_USE_STATIC_ASSERT

static_assert(lmat::is_mat_xpr<lmat::cref_matrix_ex<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::cref_matrix_ex<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::cref_matrix_ex<double> >::value, "Interface verification failed.");

static_assert(lmat::is_mat_xpr<lmat::ref_matrix_ex<double> >::value, "Interface verification failed.");
static_assert(lmat::is_mat_view<lmat::ref_matrix_ex<double> >::value, "Interface verification failed.");
static_assert(lmat::is_dense_mat<lmat::ref_matrix_ex<double> >::value, "Interface verification failed.");

#endif


MN_CASE( cref_mat_ex, constructs )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	scoped_array<double> s(ldim * n);
	const double *ps = s.ptr_begin();

	cref_matrix_ex<double, M, N> a(ps, m, n, ldim);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), ldim);
	ASSERT_EQ(a.size(), (size_t)a.nelems() );

	ASSERT_EQ(a.ptr_data(), ps);

	cref_matrix_ex<double, M, N> a2(a);

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), ldim);
	ASSERT_EQ(a2.size(), (size_t)a2.nelems() );

	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( ref_mat_ex, constructs )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	scoped_array<double> s(ldim * n);
	double *ps = s.ptr_begin();

	ref_matrix_ex<double, M, N> a(ps, m, n, ldim);

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), ldim);
	ASSERT_EQ(a.size(), (size_t)a.nelems() );

	ASSERT_EQ(a.ptr_data(), ps);

	ref_matrix_ex<double, M, N> a2(a);

	ASSERT_EQ(a2.nrows(), m);
	ASSERT_EQ(a2.ncolumns(), n);
	ASSERT_EQ(a2.nelems(), m * n);
	ASSERT_EQ(a2.lead_dim(), ldim);
	ASSERT_EQ(a2.size(), (size_t)a2.nelems() );

	ASSERT_EQ(a2.ptr_data(), ps);
}

MN_CASE( cref_mat_ex, access )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	scoped_array<double> ref(ldim * n);

	for (index_t i = 0; i < ldim * n; ++i) ref[i] = double(i + 2);
	cref_matrix_ex<double, M, N> a(ref.ptr_begin(), m, n, ldim);
	const cref_matrix_ex<double, M, N>& ac = a;

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


MN_CASE( ref_mat_ex, access )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	scoped_array<double> ref(ldim * n);

	for (index_t i = 0; i < ldim * n; ++i) ref[i] = double(i + 2);
	ref_matrix_ex<double, M, N> a(ref.ptr_begin(), m, n, ldim);
	const ref_matrix_ex<double, M, N>& ac = a;

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


MN_CASE( ref_mat_ex, assign )
{
	const index_t ldim1 = 7;
	const index_t ldim2 = 9;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	darray<double> s1(ldim1 * n);
	darray<double> s2(ldim2 * n);

	double *ps1 = s1.ptr_begin();
	double *ps2 = s2.ptr_begin();

	for (index_t i = 0; i < ldim1 * n; ++i) s1[i] = double(i + 2);
	for (index_t i = 0; i < ldim2 * n; ++i) s2[i] = double(2 * i + 3);

	darray<double> s1r(s1);

	ref_matrix_ex<double, M, N> a1(ps1, m, n, ldim1);
	ref_matrix_ex<double, M, N> a2(ps2, m, n, ldim2);

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


MN_CASE( ref_mat_ex, assign_gen )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	darray<double> ref(ldim * n);
	darray<double> s0(ldim * n);
	fill(s0, -1.0);
	darray<double> s(s0);
	double *ps = s.ptr_begin();

	// zeros

	ref_matrix_ex<double, M, N> a(ps, m, n, ldim);

	a = zeros<double>();

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim] = 0.0;
	}

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(ldim * n, ps, ref);

	// fill_value

	const double v1 = 2.5;
	a = fill_value(v1);

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim] = v1;
	}

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(ldim * n, ps, ref);

	// copy_value

	scoped_array<double> cs(m * n);

	for (index_t i = 0; i < m * n; ++i) cs[i] = double(i + 2);
	a = copy_from(cs.ptr_begin());

	ref = s0;
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
			ref[i + j * ldim] = cs[i + j * m];
	}

	ASSERT_EQ(a.nrows(), m);
	ASSERT_EQ(a.ncolumns(), n);
	ASSERT_EQ(a.nelems(), m * n);
	ASSERT_EQ(a.lead_dim(), ldim);
	ASSERT_EQ(a.ptr_data(), ps);

	ASSERT_VEC_EQ(ldim * n, ps, ref);

}



BEGIN_TPACK( cref_mat_ex_constructs )
	ADD_MN_CASE( cref_mat_ex, constructs, 0, 0 )
	ADD_MN_CASE( cref_mat_ex, constructs, 0, 1 )
	ADD_MN_CASE( cref_mat_ex, constructs, 0, 4 )
	ADD_MN_CASE( cref_mat_ex, constructs, 1, 0 )
	ADD_MN_CASE( cref_mat_ex, constructs, 1, 1 )
	ADD_MN_CASE( cref_mat_ex, constructs, 1, 4 )
	ADD_MN_CASE( cref_mat_ex, constructs, 3, 0 )
	ADD_MN_CASE( cref_mat_ex, constructs, 3, 1 )
	ADD_MN_CASE( cref_mat_ex, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_ex_constructs )
	ADD_MN_CASE( ref_mat_ex, constructs, 0, 0 )
	ADD_MN_CASE( ref_mat_ex, constructs, 0, 1 )
	ADD_MN_CASE( ref_mat_ex, constructs, 0, 4 )
	ADD_MN_CASE( ref_mat_ex, constructs, 1, 0 )
	ADD_MN_CASE( ref_mat_ex, constructs, 1, 1 )
	ADD_MN_CASE( ref_mat_ex, constructs, 1, 4 )
	ADD_MN_CASE( ref_mat_ex, constructs, 3, 0 )
	ADD_MN_CASE( ref_mat_ex, constructs, 3, 1 )
	ADD_MN_CASE( ref_mat_ex, constructs, 3, 4 )
END_TPACK

BEGIN_TPACK( cref_mat_ex_access )
	ADD_MN_CASE( cref_mat_ex, access, 0, 0 )
	ADD_MN_CASE( cref_mat_ex, access, 0, 1 )
	ADD_MN_CASE( cref_mat_ex, access, 0, 4 )
	ADD_MN_CASE( cref_mat_ex, access, 1, 0 )
	ADD_MN_CASE( cref_mat_ex, access, 1, 1 )
	ADD_MN_CASE( cref_mat_ex, access, 1, 4 )
	ADD_MN_CASE( cref_mat_ex, access, 3, 0 )
	ADD_MN_CASE( cref_mat_ex, access, 3, 1 )
	ADD_MN_CASE( cref_mat_ex, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_ex_access )
	ADD_MN_CASE( ref_mat_ex, access, 0, 0 )
	ADD_MN_CASE( ref_mat_ex, access, 0, 1 )
	ADD_MN_CASE( ref_mat_ex, access, 0, 4 )
	ADD_MN_CASE( ref_mat_ex, access, 1, 0 )
	ADD_MN_CASE( ref_mat_ex, access, 1, 1 )
	ADD_MN_CASE( ref_mat_ex, access, 1, 4 )
	ADD_MN_CASE( ref_mat_ex, access, 3, 0 )
	ADD_MN_CASE( ref_mat_ex, access, 3, 1 )
	ADD_MN_CASE( ref_mat_ex, access, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_ex_assign )
	ADD_MN_CASE( ref_mat_ex, assign, 0, 0 )
	ADD_MN_CASE( ref_mat_ex, assign, 0, 1 )
	ADD_MN_CASE( ref_mat_ex, assign, 0, 4 )
	ADD_MN_CASE( ref_mat_ex, assign, 1, 0 )
	ADD_MN_CASE( ref_mat_ex, assign, 1, 1 )
	ADD_MN_CASE( ref_mat_ex, assign, 1, 4 )
	ADD_MN_CASE( ref_mat_ex, assign, 3, 0 )
	ADD_MN_CASE( ref_mat_ex, assign, 3, 1 )
	ADD_MN_CASE( ref_mat_ex, assign, 3, 4 )
END_TPACK

BEGIN_TPACK( ref_mat_ex_assign_gen )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 0, 0 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 0, 1 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 0, 4 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 1, 0 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 1, 1 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 1, 4 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 3, 0 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 3, 1 )
	ADD_MN_CASE( ref_mat_ex, assign_gen, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( cref_mat_ex_constructs )
	ADD_TPACK( cref_mat_ex_access )

	ADD_TPACK( ref_mat_ex_constructs )
	ADD_TPACK( ref_mat_ex_access )
	ADD_TPACK( ref_mat_ex_assign )
	ADD_TPACK( ref_mat_ex_assign_gen )
END_MAIN_SUITE


