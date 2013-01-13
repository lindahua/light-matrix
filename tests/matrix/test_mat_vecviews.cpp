/**
 * @file test_mat_vecviews.cpp
 *
 * Unit testing of column/row view of a matrix
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
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


template<class Mat>
col_f64 extract_diag(const Mat& mat)
{
	const index_t m = mat.nrows();
	const index_t n = mat.ncolumns();
	const index_t len = m < n ? m : n;

	col_f64 dv(len);

	for (index_t i = 0; i < len; ++i)
	{
		dv[i] = mat(i, i);
	}

	return dv;
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

	typedef typename vecview_map<mat_t, Rgn>::const_type ccol_t;
	typedef typename vecview_map<mat_t, Rgn>::_type col_t;

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

	typedef typename vecview_map<mat_t, Rgn>::const_type crow_t;
	typedef typename vecview_map<mat_t, Rgn>::_type row_t;

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


template<template<typename T, int R, int C> class ClassT, int M, int N>
void test_diag_view()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;
	const index_t dlen = m < n ? m : n;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename diagview_map<mat_t>::const_type cdview_t;
	typedef typename diagview_map<mat_t>::_type dview_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	col_f64 dv_ref = extract_diag(a);

	cdview_t dv_c = a.diag();
	dview_t dv_u = b.diag();

	ASSERT_EQ(dv_c.nrows(), dlen);
	ASSERT_EQ(dv_c.ncolumns(), 1);
	ASSERT_EQ(dv_u.nrows(), dlen);
	ASSERT_EQ(dv_u.ncolumns(), 1);

	ASSERT_VEC_EQ(dlen, dv_c, dv_ref);
	ASSERT_VEC_EQ(dlen, dv_u, dv_ref);

}




// Columns

MN_CASE( colview_of_mat )
{
	test_col_view<ref_matrix, M, N>();
}

MN_CASE( colview_of_block )
{
	test_col_view<ref_block, M, N>();
}

MN_CASE( colview_of_grid )
{
	test_col_view<ref_grid, M, N>();
}

MN_CASE( colwhole_of_mat )
{
	test_col_range<ref_matrix, whole, M, N>();
}

MN_CASE( colwhole_of_block )
{
	test_col_range<ref_block, whole, M, N>();
}

MN_CASE( colwhole_of_grid )
{
	test_col_range<ref_grid, whole, M, N>();
}

MN_CASE( colrange_of_mat )
{
	test_col_range<ref_matrix, range, M, N>();
}

MN_CASE( colrange_of_block )
{
	test_col_range<ref_block, range, M, N>();
}

MN_CASE( colrange_of_grid )
{
	test_col_range<ref_grid, range, M, N>();
}

MN_CASE( colstep_of_mat )
{
	test_col_range<ref_matrix, step_range, M, N>();
}

MN_CASE( colstep_of_block )
{
	test_col_range<ref_block, step_range, M, N>();
}

MN_CASE( colstep_of_grid )
{
	test_col_range<ref_grid, step_range, M, N>();
}


// Rows

MN_CASE( rowview_of_mat )
{
	test_row_view<ref_matrix, M, N>();
}

MN_CASE( rowview_of_block )
{
	test_row_view<ref_block, M, N>();
}

MN_CASE( rowview_of_grid )
{
	test_row_view<ref_grid, M, N>();
}

MN_CASE( rowwhole_of_mat )
{
	test_row_range<ref_matrix, whole, M, N>();
}

MN_CASE( rowwhole_of_block )
{
	test_row_range<ref_block, whole, M, N>();
}

MN_CASE( rowwhole_of_grid )
{
	test_row_range<ref_grid, whole, M, N>();
}

MN_CASE( rowrange_of_mat )
{
	test_row_range<ref_matrix, range, M, N>();
}

MN_CASE( rowrange_of_block )
{
	test_row_range<ref_block, range, M, N>();
}

MN_CASE( rowrange_of_grid )
{
	test_row_range<ref_grid, range, M, N>();
}

MN_CASE( rowstep_of_mat )
{
	test_row_range<ref_matrix, step_range, M, N>();
}

MN_CASE( rowstep_of_block )
{
	test_row_range<ref_block, step_range, M, N>();
}

MN_CASE( rowstep_of_grid )
{
	test_row_range<ref_grid, step_range, M, N>();
}


// V-Range


// sub-vec of column

N_CASE( col_subvec_whole_of_mat )
{
	test_col_vrange<ref_matrix, whole, N>();
}

N_CASE( col_subvec_whole_of_block )
{
	test_col_vrange<ref_block, whole, N>();
}

N_CASE( col_subvec_whole_of_grid )
{
	test_col_vrange<ref_grid, whole, N>();
}

N_CASE( col_subvec_range_of_mat )
{
	test_col_vrange<ref_matrix, range, N>();
}

N_CASE( col_subvec_range_of_block )
{
	test_col_vrange<ref_block, range, N>();
}

N_CASE( col_subvec_range_of_grid )
{
	test_col_vrange<ref_grid, range, N>();
}

N_CASE( col_subvec_step_of_mat )
{
	test_col_vrange<ref_matrix, step_range, N>();
}

N_CASE( col_subvec_step_of_block )
{
	test_col_vrange<ref_block, step_range, N>();
}

N_CASE( col_subvec_step_of_grid )
{
	test_col_vrange<ref_grid, step_range, N>();
}


// sub-vec of row

N_CASE( row_subvec_whole_of_mat )
{
	test_row_vrange<ref_matrix, whole, N>();
}

N_CASE( row_subvec_whole_of_block )
{
	test_row_vrange<ref_block, whole, N>();
}

N_CASE( row_subvec_whole_of_grid )
{
	test_row_vrange<ref_grid, whole, N>();
}

N_CASE( row_subvec_range_of_mat )
{
	test_row_vrange<ref_matrix, range, N>();
}

N_CASE( row_subvec_range_of_block )
{
	test_row_vrange<ref_block, range, N>();
}

N_CASE( row_subvec_range_of_grid )
{
	test_row_vrange<ref_grid, range, N>();
}

N_CASE( row_subvec_step_of_mat )
{
	test_row_vrange<ref_matrix, step_range, N>();
}

N_CASE( row_subvec_step_of_block )
{
	test_row_vrange<ref_block, step_range, N>();
}

N_CASE( row_subvec_step_of_grid )
{
	test_row_vrange<ref_grid, step_range, N>();
}


// diag-view

MN_CASE( diagview_of_mat )
{
	test_diag_view<ref_matrix, M, N>();
}

MN_CASE( diagview_of_block )
{
	test_diag_view<ref_block, M, N>();
}

MN_CASE( diagview_of_grid )
{
	test_diag_view<ref_grid, M, N>();
}


AUTO_TPACK( colview )
{
	ADD_MN_CASE_3X3( colview_of_mat, DM, DN )
	ADD_MN_CASE_3X3( colview_of_block, DM, DN )
	ADD_MN_CASE_3X3( colview_of_grid, DM, DN )
}

AUTO_TPACK( colwhole )
{
	ADD_MN_CASE_3X3( colwhole_of_mat, DM, DN )
	ADD_MN_CASE_3X3( colwhole_of_block, DM, DN )
	ADD_MN_CASE_3X3( colwhole_of_grid, DM, DN )
}

AUTO_TPACK( colrange )
{
	ADD_MN_CASE_3X3( colrange_of_mat, DM, DN )
	ADD_MN_CASE_3X3( colrange_of_block, DM, DN )
	ADD_MN_CASE_3X3( colrange_of_grid, DM, DN )
}

AUTO_TPACK( colstep )
{
	ADD_MN_CASE_3X3( colstep_of_mat, DM, DN )
	ADD_MN_CASE_3X3( colstep_of_block, DM, DN )
	ADD_MN_CASE_3X3( colstep_of_grid, DM, DN )
}


AUTO_TPACK( rowview )
{
	ADD_MN_CASE_3X3( rowview_of_mat, DM, DN )
	ADD_MN_CASE_3X3( rowview_of_block, DM, DN )
	ADD_MN_CASE_3X3( rowview_of_grid, DM, DN )
}

AUTO_TPACK( rowwhole )
{
	ADD_MN_CASE_3X3( rowwhole_of_mat, DM, DN )
	ADD_MN_CASE_3X3( rowwhole_of_block, DM, DN )
	ADD_MN_CASE_3X3( rowwhole_of_grid, DM, DN )
}

AUTO_TPACK( rowrange )
{
	ADD_MN_CASE_3X3( rowrange_of_mat, DM, DN )
	ADD_MN_CASE_3X3( rowrange_of_block, DM, DN )
	ADD_MN_CASE_3X3( rowrange_of_grid, DM, DN )
}

AUTO_TPACK( rowstep )
{
	ADD_MN_CASE_3X3( rowstep_of_mat, DM, DN )
	ADD_MN_CASE_3X3( rowstep_of_block, DM, DN )
	ADD_MN_CASE_3X3( rowstep_of_grid, DM, DN )
}


AUTO_TPACK( diagview )
{
	ADD_MN_CASE_3X3( diagview_of_mat, DM, DN )
	ADD_MN_CASE_3X3( diagview_of_block, DM, DN )
	ADD_MN_CASE_3X3( diagview_of_grid, DM, DN )
}


AUTO_TPACK( col_subvec )
{
	ADD_N_CASE_3( col_subvec_whole_of_mat, DM )
	ADD_N_CASE_3( col_subvec_range_of_mat, DM )
	ADD_N_CASE_3( col_subvec_step_of_mat, DM )
}

AUTO_TPACK( row_subvec )
{
	ADD_N_CASE_3( row_subvec_whole_of_mat, DM )
	ADD_N_CASE_3( row_subvec_range_of_mat, DM )
	ADD_N_CASE_3( row_subvec_step_of_mat, DM )
}



