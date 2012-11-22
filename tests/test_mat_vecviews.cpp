/**
 * @file test_mat_vecviews.cpp
 *
 * Unit testing of column/row view of a matrix
 *
 * @author Dahua Lin
 */

#include "test_base.h"
#include <light_mat/matrix/matrix_classes.h>

using namespace lmat;
using namespace lmat::test;


// auxiliary functions

void fill_lin(dblock<double>& arr)
{
	for (index_t i = 0; i < arr.nelems(); ++i)
		arr[i] = double(i + 1);
}


// matrix construction facilities

const index_t DM = 7;
const index_t DN = 8;
const index_t rs = 2;
const index_t cs = 15;
const index_t LDim = 15;

template<template<typename T, int R, int C> class ClassT, int M, int N>
struct mat_maker;

template<int M, int N>
struct mat_maker<ref_matrix, M, N>
{
	typedef cref_matrix<double, M, N> cmat_t;
	typedef ref_matrix<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return m * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n);
	}
};

template<int M, int N>
struct mat_maker<ref_block, M, N>
{
	typedef cref_block<double, M, N> cmat_t;
	typedef ref_block<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return LDim * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n, cs);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n, cs);
	}
};

template<int M, int N>
struct mat_maker<ref_grid, M, N>
{
	typedef cref_grid<double, M, N> cmat_t;
	typedef ref_grid<double, M, N> mat_t;

	static index_t max_size(index_t m, index_t n)
	{
		return LDim * n;
	}

	static cmat_t get_cmat(const double *p, index_t m, index_t n)
	{
		return cmat_t(p, m, n, rs, cs);
	}

	static mat_t get_mat(double *p, index_t m, index_t n)
	{
		return mat_t(p, m, n, rs, cs);
	}
};


template<class Rgn> struct range_maker;

template<>
struct range_maker<whole>
{
	static whole get(index_t )
	{
		return whole();
	}
};

template<>
struct range_maker<range>
{
	static range get(index_t n)
	{
		return n > 2 ? range(1, n-2) : range(0, n);
	}
};

template<>
struct range_maker<step_range>
{
	static step_range get(index_t n)
	{
		if (n > 3)
		{
			index_t m = n / 2;
			return step_range(1, m, 2);
		}
		else
		{
			return step_range(0, n, 1);
		}
	}
};




// test helpers

template<class Mat, class Rgn>
col_f64 extract_col(const Mat& mat, index_t j, const Rgn& rgn)
{
	const index_t m = mat.nrows();

	const index_t len = rgn.get_num(m);
	col_f64 col(len);

	for (index_t i = 0; i < len; ++i)
	{
		col[i] = mat(rgn.get_offset(m, i), j);
	}

	return col;
}

template<class Mat, class Rgn>
row_f64 extract_row(const Mat& mat, index_t i, const Rgn& rgn)
{
	const index_t n = mat.ncolumns();

	const index_t len = rgn.get_num(n);
	row_f64 row(len);

	for (index_t j = 0; j < len; ++j)
	{
		row[j] = mat(i, rgn.get_offset(n, j));
	}

	return row;
}


template<template<typename T, int R, int C> class ClassT, int M, int N>
void test_col_view()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename colview_map<mat_t, whole>::const_type ccol_t;
	typedef typename colview_map<mat_t, whole>::non_const_type col_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	for (index_t j = 0; j < n; ++j)
	{
		col_f64 col_ref = extract_col(a, j, whole());

		ccol_t col_c = a.column(j);
		col_t col_u = b.column(j);

		ASSERT_EQ(col_c.nrows(), m);
		ASSERT_EQ(col_c.ncolumns(), 1);
		ASSERT_EQ(col_u.nrows(), m);
		ASSERT_EQ(col_u.ncolumns(), 1);

		ASSERT_VEC_EQ(m, col_c, col_ref);
		ASSERT_VEC_EQ(m, col_u, col_ref);
	}
}


template<template<typename T, int R, int C> class ClassT, class Rgn, int M, int N>
void test_col_range()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename colview_map<mat_t, Rgn>::const_type ccol_t;
	typedef typename colview_map<mat_t, Rgn>::non_const_type col_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	Rgn rgn = range_maker<Rgn>::get(m);
	const index_t len = rgn.get_num(m);

	for (index_t j = 0; j < n; ++j)
	{
		col_f64 col_ref = extract_col(a, j, rgn);

		ccol_t col_c = a(rgn, j);
		col_t col_u = b(rgn, j);

		ASSERT_EQ(col_c.nrows(), len);
		ASSERT_EQ(col_c.ncolumns(), 1);
		ASSERT_EQ(col_u.nrows(), len);
		ASSERT_EQ(col_u.ncolumns(), 1);

		ASSERT_VEC_EQ(len, col_c, col_ref);
		ASSERT_VEC_EQ(len, col_u, col_ref);
	}
}


template<template<typename T, int R, int C> class ClassT, int M, int N>
void test_row_view()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename rowview_map<mat_t, whole>::const_type crow_t;
	typedef typename rowview_map<mat_t, whole>::non_const_type row_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	for (index_t i = 0; i < m; ++i)
	{
		row_f64 row_ref = extract_row(a, i, whole());

		crow_t row_c = a.row(i);
		row_t row_u = b.row(i);

		ASSERT_EQ(row_c.nrows(), 1);
		ASSERT_EQ(row_c.ncolumns(), n);
		ASSERT_EQ(row_u.nrows(), 1);
		ASSERT_EQ(row_u.ncolumns(), n);

		ASSERT_VEC_EQ(n, row_c, row_ref);
		ASSERT_VEC_EQ(n, row_u, row_ref);
	}
}


template<template<typename T, int R, int C> class ClassT, class Rgn, int M, int N>
void test_row_range()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename rowview_map<mat_t, Rgn>::const_type crow_t;
	typedef typename rowview_map<mat_t, Rgn>::non_const_type row_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	Rgn rgn = range_maker<Rgn>::get(n);
	const index_t len = rgn.get_num(n);

	for (index_t i = 0; i < m; ++i)
	{
		row_f64 row_ref = extract_row(a, i, rgn);

		crow_t row_c = a(i, rgn);
		row_t row_u = b(i, rgn);

		ASSERT_EQ(row_c.nrows(), 1);
		ASSERT_EQ(row_c.ncolumns(), len);
		ASSERT_EQ(row_u.nrows(), 1);
		ASSERT_EQ(row_u.ncolumns(), len);

		ASSERT_VEC_EQ(len, row_c, row_ref);
		ASSERT_VEC_EQ(len, row_u, row_ref);
	}
}




MN_CASE( colview, of_mat )
{
	test_col_view<ref_matrix, M, N>();
}

MN_CASE( colview, of_block )
{
	test_col_view<ref_block, M, N>();
}

MN_CASE( colwhole, of_mat )
{
	test_col_range<ref_matrix, whole, M, N>();
}

MN_CASE( colwhole, of_block )
{
	test_col_range<ref_block, whole, M, N>();
}

MN_CASE( colrange, of_mat )
{
	test_col_range<ref_matrix, range, M, N>();
}

MN_CASE( colrange, of_block )
{
	test_col_range<ref_block, range, M, N>();
}


MN_CASE( rowview, of_mat )
{
	test_row_view<ref_matrix, M, N>();
}

MN_CASE( rowview, of_block )
{
	test_row_view<ref_block, M, N>();
}

MN_CASE( rowwhole, of_mat )
{
	test_row_range<ref_matrix, whole, M, N>();
}

MN_CASE( rowwhole, of_block )
{
	test_row_range<ref_block, whole, M, N>();
}

MN_CASE( rowrange, of_mat )
{
	test_row_range<ref_matrix, range, M, N>();
}

MN_CASE( rowrange, of_block )
{
	test_row_range<ref_block, range, M, N>();
}



BEGIN_TPACK( colview_of_mat )
	ADD_MN_CASE_3X3( colview, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colview_of_block )
	ADD_MN_CASE_3X3( colview, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colwhole_of_mat )
	ADD_MN_CASE_3X3( colwhole, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colwhole_of_block )
	ADD_MN_CASE_3X3( colwhole, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colrange_of_mat )
	ADD_MN_CASE_3X3( colrange, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colrange_of_block )
	ADD_MN_CASE_3X3( colrange, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowview_of_mat )
	ADD_MN_CASE_3X3( rowview, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowview_of_block )
	ADD_MN_CASE_3X3( rowview, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowwhole_of_mat )
	ADD_MN_CASE_3X3( rowwhole, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowwhole_of_block )
	ADD_MN_CASE_3X3( rowwhole, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowrange_of_mat )
	ADD_MN_CASE_3X3( rowrange, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowrange_of_block )
	ADD_MN_CASE_3X3( rowrange, of_block, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( colview_of_mat )
	ADD_TPACK( colview_of_block )
	ADD_TPACK( colwhole_of_mat )
	ADD_TPACK( colwhole_of_block )
	ADD_TPACK( colrange_of_mat )
	ADD_TPACK( colrange_of_block )

	ADD_TPACK( rowview_of_mat )
	ADD_TPACK( rowview_of_block )
	ADD_TPACK( rowwhole_of_mat )
	ADD_TPACK( rowwhole_of_block )
	ADD_TPACK( rowrange_of_mat )
	ADD_TPACK( rowrange_of_block )
END_MAIN_SUITE





