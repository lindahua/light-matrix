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
	typedef typename colview_map<mat_t, whole>::_type col_t;

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
	typedef typename colview_map<mat_t, Rgn>::_type col_t;

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
	typedef typename rowview_map<mat_t, whole>::type row_t;

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
	typedef typename rowview_map<mat_t, Rgn>::type row_t;

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


template<template<typename T, int R, int C> class ClassT, class Rgn, int M>
void test_col_vrange()
{
	const index_t m = M == 0 ? DM : M;

	typedef mat_maker<ClassT, M, 1> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename colview_map<mat_t, Rgn>::const_type ccol_t;
	typedef typename colview_map<mat_t, Rgn>::_type col_t;

	const index_t max_size = maker_t::max_size(m, 1);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, 1);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, 1);

	Rgn rgn = range_maker<Rgn>::get(m);
	const index_t len = rgn.get_num(m);

	col_f64 col_ref = extract_col(a, 0, rgn);

	ccol_t col_c = a(rgn);
	col_t col_u = b(rgn);

	ASSERT_EQ(col_c.nrows(), len);
	ASSERT_EQ(col_c.ncolumns(), 1);
	ASSERT_EQ(col_u.nrows(), len);
	ASSERT_EQ(col_u.ncolumns(), 1);

	ASSERT_VEC_EQ(len, col_c, col_ref);
	ASSERT_VEC_EQ(len, col_u, col_ref);
}


template<template<typename T, int R, int C> class ClassT, class Rgn, int N>
void test_row_vrange()
{
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, 1, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename rowview_map<mat_t, Rgn>::const_type crow_t;
	typedef typename rowview_map<mat_t, Rgn>::_type row_t;

	const index_t max_size = maker_t::max_size(1, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), 1, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), 1, n);

	Rgn rgn = range_maker<Rgn>::get(n);
	const index_t len = rgn.get_num(n);

	row_f64 row_ref = extract_row(a, 0, rgn);

	crow_t row_c = a(rgn);
	row_t row_u = b(rgn);

	ASSERT_EQ(row_c.nrows(), 1);
	ASSERT_EQ(row_c.ncolumns(), len);
	ASSERT_EQ(row_u.nrows(), 1);
	ASSERT_EQ(row_u.ncolumns(), len);

	ASSERT_VEC_EQ(len, row_c, row_ref);
	ASSERT_VEC_EQ(len, row_u, row_ref);
}




// Columns

MN_CASE( colview, of_mat )
{
	test_col_view<ref_matrix, M, N>();
}

MN_CASE( colview, of_block )
{
	test_col_view<ref_block, M, N>();
}

MN_CASE( colview, of_grid )
{
	test_col_view<ref_grid, M, N>();
}

MN_CASE( colwhole, of_mat )
{
	test_col_range<ref_matrix, whole, M, N>();
}

MN_CASE( colwhole, of_block )
{
	test_col_range<ref_block, whole, M, N>();
}

MN_CASE( colwhole, of_grid )
{
	test_col_range<ref_grid, whole, M, N>();
}

MN_CASE( colrange, of_mat )
{
	test_col_range<ref_matrix, range, M, N>();
}

MN_CASE( colrange, of_block )
{
	test_col_range<ref_block, range, M, N>();
}

MN_CASE( colrange, of_grid )
{
	test_col_range<ref_grid, range, M, N>();
}

MN_CASE( colstep, of_mat )
{
	test_col_range<ref_matrix, step_range, M, N>();
}

MN_CASE( colstep, of_block )
{
	test_col_range<ref_block, step_range, M, N>();
}

MN_CASE( colstep, of_grid )
{
	test_col_range<ref_grid, step_range, M, N>();
}


// Rows

MN_CASE( rowview, of_mat )
{
	test_row_view<ref_matrix, M, N>();
}

MN_CASE( rowview, of_block )
{
	test_row_view<ref_block, M, N>();
}

MN_CASE( rowview, of_grid )
{
	test_row_view<ref_grid, M, N>();
}

MN_CASE( rowwhole, of_mat )
{
	test_row_range<ref_matrix, whole, M, N>();
}

MN_CASE( rowwhole, of_block )
{
	test_row_range<ref_block, whole, M, N>();
}

MN_CASE( rowwhole, of_grid )
{
	test_row_range<ref_grid, whole, M, N>();
}

MN_CASE( rowrange, of_mat )
{
	test_row_range<ref_matrix, range, M, N>();
}

MN_CASE( rowrange, of_block )
{
	test_row_range<ref_block, range, M, N>();
}

MN_CASE( rowrange, of_grid )
{
	test_row_range<ref_grid, range, M, N>();
}

MN_CASE( rowstep, of_mat )
{
	test_row_range<ref_matrix, step_range, M, N>();
}

MN_CASE( rowstep, of_block )
{
	test_row_range<ref_block, step_range, M, N>();
}

MN_CASE( rowstep, of_grid )
{
	test_row_range<ref_grid, step_range, M, N>();
}


// V-Range


// sub-vec of column

N_CASE( col_subvec, whole_of_mat )
{
	test_col_vrange<ref_matrix, whole, N>();
}

N_CASE( col_subvec, whole_of_block )
{
	test_col_vrange<ref_block, whole, N>();
}

N_CASE( col_subvec, whole_of_grid )
{
	test_col_vrange<ref_grid, whole, N>();
}

N_CASE( col_subvec, range_of_mat )
{
	test_col_vrange<ref_matrix, range, N>();
}

N_CASE( col_subvec, range_of_block )
{
	test_col_vrange<ref_block, range, N>();
}

N_CASE( col_subvec, range_of_grid )
{
	test_col_vrange<ref_grid, range, N>();
}

N_CASE( col_subvec, step_of_mat )
{
	test_col_vrange<ref_matrix, step_range, N>();
}

N_CASE( col_subvec, step_of_block )
{
	test_col_vrange<ref_block, step_range, N>();
}

N_CASE( col_subvec, step_of_grid )
{
	test_col_vrange<ref_grid, step_range, N>();
}


// sub-vec of row

N_CASE( row_subvec, whole_of_mat )
{
	test_row_vrange<ref_matrix, whole, N>();
}

N_CASE( row_subvec, whole_of_block )
{
	test_row_vrange<ref_block, whole, N>();
}

N_CASE( row_subvec, whole_of_grid )
{
	test_row_vrange<ref_grid, whole, N>();
}

N_CASE( row_subvec, range_of_mat )
{
	test_row_vrange<ref_matrix, range, N>();
}

N_CASE( row_subvec, range_of_block )
{
	test_row_vrange<ref_block, range, N>();
}

N_CASE( row_subvec, range_of_grid )
{
	test_row_vrange<ref_grid, range, N>();
}

N_CASE( row_subvec, step_of_mat )
{
	test_row_vrange<ref_matrix, step_range, N>();
}

N_CASE( row_subvec, step_of_block )
{
	test_row_vrange<ref_block, step_range, N>();
}

N_CASE( row_subvec, step_of_grid )
{
	test_row_vrange<ref_grid, step_range, N>();
}




BEGIN_TPACK( colview_of_mat )
	ADD_MN_CASE_3X3( colview, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colview_of_block )
	ADD_MN_CASE_3X3( colview, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colview_of_grid )
	ADD_MN_CASE_3X3( colview, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( colwhole_of_mat )
	ADD_MN_CASE_3X3( colwhole, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colwhole_of_block )
	ADD_MN_CASE_3X3( colwhole, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colwhole_of_grid )
	ADD_MN_CASE_3X3( colwhole, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( colrange_of_mat )
	ADD_MN_CASE_3X3( colrange, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colrange_of_block )
	ADD_MN_CASE_3X3( colrange, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colrange_of_grid )
	ADD_MN_CASE_3X3( colrange, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( colstep_of_mat )
	ADD_MN_CASE_3X3( colstep, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( colstep_of_block )
	ADD_MN_CASE_3X3( colstep, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( colstep_of_grid )
	ADD_MN_CASE_3X3( colstep, of_grid, DM, DN )
END_TPACK


BEGIN_TPACK( rowview_of_mat )
	ADD_MN_CASE_3X3( rowview, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowview_of_block )
	ADD_MN_CASE_3X3( rowview, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowview_of_grid )
	ADD_MN_CASE_3X3( rowview, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( rowwhole_of_mat )
	ADD_MN_CASE_3X3( rowwhole, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowwhole_of_block )
	ADD_MN_CASE_3X3( rowwhole, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowwhole_of_grid )
	ADD_MN_CASE_3X3( rowwhole, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( rowrange_of_mat )
	ADD_MN_CASE_3X3( rowrange, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowrange_of_block )
	ADD_MN_CASE_3X3( rowrange, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowrange_of_grid )
	ADD_MN_CASE_3X3( rowrange, of_grid, DM, DN )
END_TPACK

BEGIN_TPACK( rowstep_of_mat )
	ADD_MN_CASE_3X3( rowstep, of_mat, DM, DN )
END_TPACK

BEGIN_TPACK( rowstep_of_block )
	ADD_MN_CASE_3X3( rowstep, of_block, DM, DN )
END_TPACK

BEGIN_TPACK( rowstep_of_grid )
	ADD_MN_CASE_3X3( rowstep, of_grid, DM, DN )
END_TPACK


BEGIN_TPACK( col_subvec )
	ADD_N_CASE( col_subvec, whole_of_mat, 0 )
	ADD_N_CASE( col_subvec, whole_of_block, 1 )
	ADD_N_CASE( col_subvec, whole_of_grid, DM )

	ADD_N_CASE( col_subvec, range_of_mat, 0 )
	ADD_N_CASE( col_subvec, range_of_block, 1 )
	ADD_N_CASE( col_subvec, range_of_grid, DM )

	ADD_N_CASE( col_subvec, step_of_mat, 0 )
	ADD_N_CASE( col_subvec, step_of_block, 1 )
	ADD_N_CASE( col_subvec, step_of_grid, DM )
END_TPACK

BEGIN_TPACK( row_subvec )
	ADD_N_CASE( row_subvec, whole_of_mat, 0 )
	ADD_N_CASE( row_subvec, whole_of_grid, DM )

	ADD_N_CASE( row_subvec, range_of_mat, 0 )
	ADD_N_CASE( row_subvec, range_of_grid, DM )

	ADD_N_CASE( row_subvec, step_of_mat, 0 )
	ADD_N_CASE( row_subvec, step_of_grid, DM )
END_TPACK



BEGIN_MAIN_SUITE
	// col views

	ADD_TPACK( colview_of_mat )
	ADD_TPACK( colview_of_block )
	ADD_TPACK( colview_of_grid )

	ADD_TPACK( colwhole_of_mat )
	ADD_TPACK( colwhole_of_block )
	ADD_TPACK( colwhole_of_grid )

	ADD_TPACK( colrange_of_mat )
	ADD_TPACK( colrange_of_block )
	ADD_TPACK( colrange_of_grid )

	ADD_TPACK( colstep_of_mat )
	ADD_TPACK( colstep_of_block )
	ADD_TPACK( colstep_of_grid )

	// row views

	ADD_TPACK( rowview_of_mat )
	ADD_TPACK( rowview_of_block )
	ADD_TPACK( rowview_of_grid )

	ADD_TPACK( rowwhole_of_mat )
	ADD_TPACK( rowwhole_of_block )
	ADD_TPACK( rowwhole_of_grid )

	ADD_TPACK( rowrange_of_mat )
	ADD_TPACK( rowrange_of_block )
	ADD_TPACK( rowrange_of_grid )

	ADD_TPACK( rowstep_of_mat )
	ADD_TPACK( rowstep_of_block )
	ADD_TPACK( rowstep_of_grid )

	// sub vectors

	ADD_TPACK( col_subvec )
	ADD_TPACK( row_subvec )
END_MAIN_SUITE





