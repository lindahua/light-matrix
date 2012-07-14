/**
 * @file test_dense_veval.cpp
 *
 * Unit testing of vector evaluators for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matrix/generic_matrix_eval.h>

using namespace lmat;
using namespace lmat::test;

#ifdef LMAT_USE_STATIC_ASSERT
static_assert(is_linear_vector_evaluator<continuous_linear_evaluator<double>, double>::value,
		"Evaluator interface check failed");
static_assert(is_linear_vector_evaluator<cached_linear_evaluator<double>, double>::value,
		"Evaluator interface check failed");
static_assert(is_linear_vector_evaluator<const_linear_evaluator<double>, double>::value,
		"Evaluator interface check failed");

static_assert(is_percol_vector_evaluator<dense_percol_evaluator<double>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<cached_percol_evaluator<double>, double>::value,
		"Evaluator interface check failed");
static_assert(is_percol_vector_evaluator<const_percol_evaluator<double>, double>::value,
		"Evaluator interface check failed");
#endif


MN_CASE( linear_veval, continu_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef dense_matrix<double, M, N> mat;

	typedef continuous_linear_evaluator<double> veval;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert(is_same<typename linear_eval<mat>::evaluator_type, veval>::value,
			"Evaluator type verification failed");
#endif

	mat a(m, n);
	veval ve(a);

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ve.get_value(i), a[i]);
	}
}


MN_CASE( linear_veval, cached_linear )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;

	typedef continuous_linear_evaluator<double> veval_ex1;
	typedef cached_linear_evaluator<double> veval_ex;

#ifdef LMAT_USE_STATIC_ASSERT

	static_assert(
			(N == 1 && is_same<typename linear_eval<mat_ex>::evaluator_type, veval_ex1>::value) ||
			(N != 1  && is_same<typename linear_eval<mat_ex>::evaluator_type, veval_ex>::value),
			"Evaluator type verification failed");
#endif

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);

	veval_ex ve(a);

	dense_matrix<double, M, N> ar(a);

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ve.get_value(i), ar[i]);
	}
}


MN_CASE( linear_veval, const_linear )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef const_matrix<double, M, N> mat;
	typedef const_linear_evaluator<double> veval;

#ifdef LMAT_USE_STATIC_ASSERT
	static_assert(is_same<typename linear_eval<mat>::evaluator_type, veval>::value,
			"Evaluator type verification failed");
#endif

	const double v = 12.5;
	mat a(m, n, v);
	veval ve(a);

	for (index_t i = 0; i < m * n; ++i)
	{
		ASSERT_EQ(ve.get_value(i), v);
	}
}


MN_CASE( percol_veval, dense_percol )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;
	typedef dense_percol_evaluator<double> veval_ex;

#ifdef LMAT_USE_STATIC_ASSERT

	static_assert(
			is_same<typename percol_eval<mat_ex>::evaluator_type, veval_ex>::value,
			"Evaluator type verification failed");
#endif

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);

	veval_ex ve(a);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ve.get_value(i), a(i, j));
		}

		ve.next_column();
	}
}


MN_CASE( percol_veval, cached_percol )
{
	const index_t ldim = 7;
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef ref_matrix_ex<double, M, N> mat_ex;
	typedef cached_percol_evaluator<double> veval_c;

	scoped_array<double> s(ldim * n);
	mat_ex a(s.ptr_begin(), m, n, ldim);

	veval_c ve(a);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ve.get_value(i), a(i, j));
		}

		ve.next_column();
	}
}


MN_CASE( percol_veval, const_percol )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	typedef const_matrix<double, M, N> mat;
	typedef const_percol_evaluator<double> veval;

#ifdef LMAT_USE_STATIC_ASSERT

	static_assert(
			is_same<typename percol_eval<mat>::evaluator_type, veval>::value,
			"Evaluator type verification failed");
#endif

	const double val = 12.5;
	mat a(m, n, val);
	veval ve(a);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			ASSERT_EQ(ve.get_value(i), val);
		}

		ve.next_column();
	}
}


MN_CASE( eval_by_scalars, cont_to_cont )
{
	const index_t m = M == 0 ? 5 : M;
	const index_t n = N == 0 ? 6 : N;

	dense_matrix<double, M, N> a(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = double(i+2);

	dense_matrix<double, M, N> b(m, n, zeros<double>());

	evaluate_by_scalars(a, b);

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

	evaluate_by_scalars(a, b);

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

	evaluate_by_scalars(a, b);

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

	evaluate_by_scalars(a, b);

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



