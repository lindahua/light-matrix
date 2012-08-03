/**
 * @file test_dense_veval.cpp
 *
 * Unit testing of vector evaluators for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_vector_eval.h>

using namespace lmat;
using namespace lmat::test;


MN_CASE( linear_veval, continu_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, N> mat;

	typedef continuous_linear_veval_scheme<double> sch_t;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert(is_same<
			typename vector_eval<mat, as_linear_vec, scalar_kernel_t>::scheme_type,
			sch_t>::value,
			"Scheme type verification failed");
#endif

	mat a(m, n);
	sch_t sch(a);
	const IVecEvalLinearScheme<sch_t, double>& sch_r = sch;

	veval_memacc_kernel<double> ker = sch_r.kernel();
	const double *st = sch_r.vec_state();

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ker.get_value(i, st), a[i]);
	}
}


MN_CASE( linear_veval, cached_linear )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;

	typedef typename
			if_c<N == 1,
			continuous_linear_veval_scheme<double>,
			cached_linear_veval_scheme<double, M, N> >::type sch_t;

#ifdef LMAT_USE_STATIC_ASSERT

	static_assert(is_same<
			typename vector_eval<mat_ex, as_linear_vec, scalar_kernel_t>::scheme_type,
			sch_t>::value,
			"Scheme type verification failed");
#endif

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);

	sch_t sch(a);
	const IVecEvalLinearScheme<sch_t, double>& sch_r = sch;

	veval_memacc_kernel<double> ker = sch_r.kernel();
	const double *st = sch_r.vec_state();

	dense_matrix<double, M, N> ar(a);

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ker.get_value(i, st), ar[i]);
	}
}


MN_CASE( linear_veval, const_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef const_matrix<double, M, N> mat;
	typedef const_linear_veval_scheme<double> sch_t;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert(is_same<
			typename vector_eval<mat, as_linear_vec, scalar_kernel_t>::scheme_type,
			sch_t>::value,
			"Scheme type verification failed");
#endif

	const double v = 12.5;
	mat a(m, n, v);
	sch_t sch(a);
	const IVecEvalLinearScheme<sch_t, double>& sch_r = sch;

	veval_const_kernel<double> ker = sch_r.kernel();
	double st = sch_r.vec_state();

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ker.get_value(i, st), v);
	}
}


MN_CASE( percol_veval, dense_percol )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;
	typedef dense_percol_veval_scheme<double> sch_t;

#ifdef LMAT_USE_STATIC_ASSERT

	static_assert(is_same<
			typename vector_eval<mat_ex, per_column, scalar_kernel_t>::scheme_type,
			sch_t>::value,
			"Scheme type verification failed");
#endif

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);
	sch_t sch(a);
	const IVecEvalPerColScheme<sch_t, double>& sch_r = sch;

	veval_memacc_kernel<double> ker = sch_r.kernel();

	for (index_t j = 0; j < n; ++j)
	{
		const double *st = sch_r.col_state(j);
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ker.get_value(i, st), a(i, j));
		}
	}
}


MN_CASE( percol_veval, cached_percol )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;
	typedef cached_percol_veval_scheme<double, M, N> sch_t;

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);
	sch_t sch(a);
	const IVecEvalPerColScheme<sch_t, double>& sch_r = sch;

	veval_memacc_kernel<double> ker = sch_r.kernel();

	for (index_t j = 0; j < n; ++j)
	{
		const double* s = sch_r.col_state(j);

		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ker.get_value(i, s), a(i, j));
		}
	}
}


MN_CASE( percol_veval, const_percol )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef const_matrix<double, M, N> mat;
	typedef const_percol_veval_scheme<double> sch_t;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert(is_same<
			typename vector_eval<mat, per_column, scalar_kernel_t>::scheme_type,
			sch_t>::value,
			"Scheme type verification failed");
#endif

	const double val = 12.5;
	mat a(m, n, val);
	sch_t sch(a);
	const IVecEvalPerColScheme<sch_t, double>& sch_r = sch;

	veval_const_kernel<double> ker = sch_r.kernel();

	for (index_t j = 0; j < n; ++j)
	{
		double st = sch_r.col_state(j);
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ker.get_value(i, st), val);
		}
	}
}



MN_CASE( eval_by_scalars, cont_to_cont )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> a(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = double(i+2);

	dense_matrix<double, M, N> b(m, n, zeros<double>());

	linear_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);

	fill(b, 0.0);

	percol_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);
}

MN_CASE( eval_by_scalars, cont_to_ext )
{
	const index_t ldim_b = 8;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> a(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = double(i+2);

	scoped_array<double> sb(ldim_b * n);
	fill(sb, 0.0);
	ref_matrix_ex<double, M, N> b(sb.ptr_begin(), m, n, ldim_b);

	percol_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);
}

MN_CASE( eval_by_scalars, ext_to_cont )
{
	const index_t ldim_a = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<double> sa(ldim_a * n);
	for (index_t i = 0; i < ldim_a * n; ++i) sa[i] = double(i+2);

	ref_matrix_ex<double, M, N> a(sa.ptr_begin(), m, n, ldim_a);
	dense_matrix<double, M, N> b(m, n, fill_value(0.0));

	linear_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);

	fill(b, 0.0);

	percol_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);
}

MN_CASE( eval_by_scalars, ext_to_ext )
{
	const index_t ldim_a = 7;
	const index_t ldim_b = 8;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	scoped_array<double> sa(ldim_a * n);
	for (index_t i = 0; i < ldim_a * n; ++i) sa[i] = double(i+2);

	scoped_array<double> sb(ldim_b * n);
	fill(sb, 0.0);

	ref_matrix_ex<double, M, N> a(sa.ptr_begin(), m, n, ldim_a);
	ref_matrix_ex<double, M, N> b(sb.ptr_begin(), m, n, ldim_b);

	percol_by_scalars_evaluate(a, b);
	ASSERT_MAT_EQ(m, n, a, b);
}



BEGIN_TPACK( continu_linear_eval )
	ADD_MN_CASE_3X3( linear_veval, continu_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( cached_linear_eval )
	ADD_MN_CASE_3X3( linear_veval, cached_linear, 5, 6 )
END_TPACK

BEGIN_TPACK( const_linear_eval )
	ADD_MN_CASE_3X3( linear_veval, const_linear, 5, 6 )
END_TPACK


BEGIN_TPACK( dense_percol_eval )
	ADD_MN_CASE_3X3( percol_veval, dense_percol, 5, 6 )
END_TPACK

BEGIN_TPACK( cached_percol_eval )
	ADD_MN_CASE_3X3( percol_veval, cached_percol, 5, 6 )
END_TPACK

BEGIN_TPACK( const_percol_eval )
	ADD_MN_CASE_3X3( percol_veval, const_percol, 5, 6 )
END_TPACK


BEGIN_TPACK( eval_by_scalars_c2c )
	ADD_MN_CASE_3X3( eval_by_scalars, cont_to_cont, 5, 6 )
END_TPACK

BEGIN_TPACK( eval_by_scalars_c2e )
	ADD_MN_CASE_3X3( eval_by_scalars, cont_to_ext, 5, 6 )
END_TPACK

BEGIN_TPACK( eval_by_scalars_e2c )
	ADD_MN_CASE_3X3( eval_by_scalars, ext_to_cont, 5, 6 )
END_TPACK

BEGIN_TPACK( eval_by_scalars_e2e )
	ADD_MN_CASE_3X3( eval_by_scalars, ext_to_ext, 5, 6 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( continu_linear_eval )
	ADD_TPACK( cached_linear_eval )
	ADD_TPACK( const_linear_eval )

	ADD_TPACK( dense_percol_eval )
	ADD_TPACK( cached_percol_eval )
	ADD_TPACK( const_percol_eval )

	ADD_TPACK( eval_by_scalars_c2c )
	ADD_TPACK( eval_by_scalars_c2e )
	ADD_TPACK( eval_by_scalars_e2c )
	ADD_TPACK( eval_by_scalars_e2e )
END_MAIN_SUITE



