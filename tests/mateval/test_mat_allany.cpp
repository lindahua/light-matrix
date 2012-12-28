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


const index_t max_len = 28;


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




T_CASE( full_reduc, all_true )
{
	dense_col<T> a(max_len);
	dense_col<T> b(max_len);

	for (index_t i = 0; i < max_len; ++i) a[i] = T(i+1);

	for (index_t k = 0; k <= max_len; ++k)
	{
		auto ak = a(range(0, k));
		auto bk = b(range(0, k));

		// all-eq
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i];
		ASSERT_EQ( all(ak == bk, true),  my_all_eq(ak, bk) );

		// all-ne
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + 1;
		ASSERT_EQ( all(ak == bk, true),  my_all_eq(ak, bk) );

		// half-half
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + T(i % 2);
		ASSERT_EQ( all(ak == bk, true),  my_all_eq(ak, bk) );
	}
}

T_CASE( full_reduc, all_false )
{
	dense_col<T> a(max_len);
	dense_col<T> b(max_len);

	for (index_t i = 0; i < max_len; ++i) a[i] = T(i+1);

	for (index_t k = 0; k <= max_len; ++k)
	{
		auto ak = a(range(0, k));
		auto bk = b(range(0, k));

		// all-eq
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i];
		ASSERT_EQ( all(ak == bk, false),  my_all_ne(ak, bk) );

		// all-ne
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + 1;
		ASSERT_EQ( all(ak == bk, false),  my_all_ne(ak, bk) );

		// half-half
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + T(i % 2);
		ASSERT_EQ( all(ak == bk, false),  my_all_ne(ak, bk) );
	}
}

T_CASE( full_reduc, any_true )
{
	dense_col<T> a(max_len);
	dense_col<T> b(max_len);

	for (index_t i = 0; i < max_len; ++i) a[i] = T(i+1);

	for (index_t k = 0; k <= max_len; ++k)
	{
		auto ak = a(range(0, k));
		auto bk = b(range(0, k));

		// all-eq
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i];
		ASSERT_EQ( any(ak == bk, true),  my_any_eq(ak, bk) );

		// all-ne
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + 1;
		ASSERT_EQ( any(ak == bk, true),  my_any_eq(ak, bk) );

		// half-half
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + T(i % 2);
		ASSERT_EQ( any(ak == bk, true),  my_any_eq(ak, bk) );
	}
}

T_CASE( full_reduc, any_false )
{
	dense_col<T> a(max_len);
	dense_col<T> b(max_len);

	for (index_t i = 0; i < max_len; ++i) a[i] = T(i+1);

	for (index_t k = 0; k <= max_len; ++k)
	{
		auto ak = a(range(0, k));
		auto bk = b(range(0, k));

		// all-eq
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i];
		ASSERT_EQ( any(ak == bk, false),  my_any_ne(ak, bk) );

		// all-ne
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + 1;
		ASSERT_EQ( any(ak == bk, false),  my_any_ne(ak, bk));

		// half-half
		for (index_t i = 0; i < k; ++i) bk[i] = ak[i] + T(i % 2);
		ASSERT_EQ( any(ak == bk, false),  my_any_ne(ak, bk) );
	}
}



BEGIN_TPACK( mat_full_allany )
	ADD_T_CASE( full_reduc, all_true,  double )
	ADD_T_CASE( full_reduc, all_true,  int )
	ADD_T_CASE( full_reduc, all_false, double )
	ADD_T_CASE( full_reduc, all_false, int )
	ADD_T_CASE( full_reduc, any_true,  double )
	ADD_T_CASE( full_reduc, any_true,  int )
	ADD_T_CASE( full_reduc, any_false, double )
	ADD_T_CASE( full_reduc, any_false, int )
END_TPACK

BEGIN_MAIN_SUITE
	ADD_TPACK( mat_full_allany )
END_MAIN_SUITE


