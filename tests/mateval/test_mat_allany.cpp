/**
 * @file test_mat_allany.cpp
 *
 * @brief Unit testing of any all reduction
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_allany.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 5;
const index_t DN = 8;

template<class A, class B>
inline bool my_all_eq(const A& a, const B& b)
{
	index_t m = common_nrows(a, b);
	index_t n = common_ncols(a, b);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (a(i, j) != b(i, j)) return false;
		}
	}

	return true;
}

template<class A, class B>
inline bool my_all_ne(const A& a, const B& b)
{
	index_t m = common_nrows(a, b);
	index_t n = common_ncols(a, b);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (a(i, j) == b(i, j)) return false;
		}
	}

	return true;
}


template<class A, class B>
inline bool my_any_eq(const A& a, const B& b)
{
	index_t m = common_nrows(a, b);
	index_t n = common_ncols(a, b);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (a(i, j) == b(i, j)) return true;
		}
	}

	return false;
}


template<class A, class B>
inline bool my_any_ne(const A& a, const B& b)
{
	index_t m = common_nrows(a, b);
	index_t n = common_ncols(a, b);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (a(i, j) != b(i, j)) return true;
		}
	}

	return false;
}


TMN_CASE( full_reduc, all_true )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	// all-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	ASSERT_EQ( all(a == b, true),  my_all_eq(a, b) );
	ASSERT_EQ( all(to_bool(a == b), true),  my_all_eq(a, b) );

	// all-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	ASSERT_EQ( all(a == b, true),  my_all_eq(a, b) );
	ASSERT_EQ( all(to_bool(a == b), true),  my_all_eq(a, b) );

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	ASSERT_EQ( all(a == b, true),  my_all_eq(a, b) );
	ASSERT_EQ( all(to_bool(a == b), true),  my_all_eq(a, b) );
}

TMN_CASE( full_reduc, all_false )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	// all-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	ASSERT_EQ( all(a == b, false),  my_all_ne(a, b) );
	ASSERT_EQ( all(to_bool(a == b), false),  my_all_ne(a, b) );

	// all-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	ASSERT_EQ( all(a == b, false),  my_all_ne(a, b) );
	ASSERT_EQ( all(to_bool(a == b), false),  my_all_ne(a, b) );

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	ASSERT_EQ( all(a == b, false),  my_all_ne(a, b) );
	ASSERT_EQ( all(to_bool(a == b), false),  my_all_ne(a, b) );
}

TMN_CASE( full_reduc, any_true )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	// any-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	ASSERT_EQ( any(a == b, true),  my_any_eq(a, b) );
	ASSERT_EQ( any(to_bool(a == b), true),  my_any_eq(a, b) );

	// any-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	ASSERT_EQ( any(a == b, true),  my_any_eq(a, b) );
	ASSERT_EQ( any(to_bool(a == b), true),  my_any_eq(a, b) );

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	ASSERT_EQ( any(a == b, true),  my_any_eq(a, b) );
	ASSERT_EQ( any(to_bool(a == b), true),  my_any_eq(a, b) );
}

TMN_CASE( full_reduc, any_false )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);

	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	// any-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	ASSERT_EQ( any(a == b, false),  my_any_ne(a, b) );
	ASSERT_EQ( any(to_bool(a == b), false),  my_any_ne(a, b) );

	// any-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	ASSERT_EQ( any(a == b, false),  my_any_ne(a, b) );
	ASSERT_EQ( any(to_bool(a == b), false),  my_any_ne(a, b) );

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	ASSERT_EQ( any(a == b, false),  my_any_ne(a, b) );
	ASSERT_EQ( any(to_bool(a == b), false),  my_any_ne(a, b) );
}


TMN_CASE( colwise_reduc, all_true )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	dense_row<bool, N> d(n);
	dense_row<bool, N> r(n);

	// all-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	for (index_t j = 0; j < n; ++j) r[j] = my_all_eq(a.column(j), b.column(j));
	colwise_all(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);

	// all-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	for (index_t j = 0; j < n; ++j) r[j] = my_all_eq(a.column(j), b.column(j));
	colwise_all(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	for (index_t j = 0; j < n; ++j) r[j] = my_all_eq(a.column(j), b.column(j));
	colwise_all(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);
}


TMN_CASE( colwise_reduc, all_false )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	dense_row<bool, N> d(n);
	dense_row<bool, N> r(n);

	// all-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	for (index_t j = 0; j < n; ++j) r[j] = my_all_ne(a.column(j), b.column(j));
	colwise_all(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

	// all-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	for (index_t j = 0; j < n; ++j) r[j] = my_all_ne(a.column(j), b.column(j));
	colwise_all(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	for (index_t j = 0; j < n; ++j) r[j] = my_all_ne(a.column(j), b.column(j));
	colwise_all(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_all(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

}


TMN_CASE( colwise_reduc, any_true )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	dense_row<bool, N> d(n);
	dense_row<bool, N> r(n);

	// any-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	for (index_t j = 0; j < n; ++j) r[j] = my_any_eq(a.column(j), b.column(j));
	colwise_any(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);

	// any-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	for (index_t j = 0; j < n; ++j) r[j] = my_any_eq(a.column(j), b.column(j));
	colwise_any(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	for (index_t j = 0; j < n; ++j) r[j] = my_any_eq(a.column(j), b.column(j));
	colwise_any(a == b, d, true);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, true);
	ASSERT_VEC_EQ(n, d, r);
}


TMN_CASE( colwise_reduc, any_false )
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	dense_matrix<T, M, N> a(m, n);
	dense_matrix<T, M, N> b(m, n);
	for (index_t i = 0; i < m * n; ++i) a[i] = T(i+1);

	dense_row<bool, N> d(n);
	dense_row<bool, N> r(n);

	// any-eq
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i];
	for (index_t j = 0; j < n; ++j) r[j] = my_any_ne(a.column(j), b.column(j));
	colwise_any(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

	// any-ne
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + 1;
	for (index_t j = 0; j < n; ++j) r[j] = my_any_ne(a.column(j), b.column(j));
	colwise_any(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

	// half-half
	for (index_t i = 0; i < m * n; ++i) b[i] = a[i] + T(i % 2);
	for (index_t j = 0; j < n; ++j) r[j] = my_any_ne(a.column(j), b.column(j));
	colwise_any(a == b, d, false);
	ASSERT_VEC_EQ(n, d, r);
	colwise_any(to_bool(a == b), d, false);
	ASSERT_VEC_EQ(n, d, r);

}



BEGIN_TPACK( mat_full_alltrue )
	ADD_TMN_CASE_3X3( full_reduc, all_true, double, DM, DN )
	ADD_TMN_CASE_3X3( full_reduc, all_true, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_full_allfalse )
	ADD_TMN_CASE_3X3( full_reduc, all_false, double, DM, DN )
	ADD_TMN_CASE_3X3( full_reduc, all_false, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_full_anytrue )
	ADD_TMN_CASE_3X3( full_reduc, any_true, double, DM, DN )
	ADD_TMN_CASE_3X3( full_reduc, any_true, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_full_anyfalse )
	ADD_TMN_CASE_3X3( full_reduc, any_false, double, DM, DN )
	ADD_TMN_CASE_3X3( full_reduc, any_false, int, DM, DN )
END_TPACK


BEGIN_TPACK( mat_colwise_alltrue )
	ADD_TMN_CASE_3X3( colwise_reduc, all_true, double, DM, DN )
	ADD_TMN_CASE_3X3( colwise_reduc, all_true, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_colwise_allfalse )
	ADD_TMN_CASE_3X3( colwise_reduc, all_false, double, DM, DN )
	ADD_TMN_CASE_3X3( colwise_reduc, all_false, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_colwise_anytrue )
	ADD_TMN_CASE_3X3( colwise_reduc, any_true, double, DM, DN )
	ADD_TMN_CASE_3X3( colwise_reduc, any_true, int, DM, DN )
END_TPACK

BEGIN_TPACK( mat_colwise_anyfalse )
	ADD_TMN_CASE_3X3( colwise_reduc, any_false, double, DM, DN )
	ADD_TMN_CASE_3X3( colwise_reduc, any_false, int, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_full_alltrue )
	ADD_TPACK( mat_full_allfalse )
	ADD_TPACK( mat_full_anytrue )
	ADD_TPACK( mat_full_anyfalse )

	ADD_TPACK( mat_colwise_alltrue )
	ADD_TPACK( mat_colwise_allfalse )
	ADD_TPACK( mat_colwise_anytrue )
	ADD_TPACK( mat_colwise_anyfalse )
END_MAIN_SUITE


