/**
 * @file test_mat_matviews.cpp
 *
 * @brief Unit testing of matrix subviews
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

template<class Mat, class RowRgn, class ColRgn>
mat_f64 extract_submat(const Mat& mat, const RowRgn& r, const ColRgn& c)
{
	const index_t m0 = mat.nrows();
	const index_t n0 = mat.ncolumns();

	const index_t sm = r.get_num(m0);
	const index_t sn = c.get_num(n0);

	mat_f64 sub(sm, sn);

	for (index_t j = 0; j < sn; ++j)
	{
		for (index_t i = 0; i < sm; ++i)
		{
			sub(i, j) = mat(r.get_offset(m0, i), c.get_offset(n0, j));
		}
	}

	return sub;
}


template<template<typename T, int R, int C> class ClassT,
	class RowRgn, class ColRgn, int M, int N>
void test_mat_range()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	typedef mat_maker<ClassT, M, N> maker_t;
	typedef typename maker_t::cmat_t cmat_t;
	typedef typename maker_t::mat_t mat_t;

	typedef typename matview_map<mat_t, RowRgn, ColRgn>::const_type csub_t;
	typedef typename matview_map<mat_t, RowRgn, ColRgn>::_type sub_t;

	const index_t max_size = maker_t::max_size(m, n);

	dblock<double> s(max_size, zero());
	fill_lin(s);

	cmat_t a = maker_t::get_cmat(s.ptr_data(), m, n);
	mat_t b = maker_t::get_mat(s.ptr_data(), m, n);

	RowRgn r = range_maker<RowRgn>::get(m);
	ColRgn c = range_maker<ColRgn>::get(n);

	const index_t sm = r.get_num(m);
	const index_t sn = c.get_num(n);

	mat_f64 sub_ref = extract_submat(a, r, c);

	csub_t sub_c = a(r, c);
	sub_t sub_u = b(r, c);

	ASSERT_EQ(sub_c.nrows(), sm);
	ASSERT_EQ(sub_c.ncolumns(), sn);
	ASSERT_EQ(sub_u.nrows(), sm);
	ASSERT_EQ(sub_u.ncolumns(), sn);

	ASSERT_MAT_EQ(sm, sn, sub_c, sub_ref);
	ASSERT_MAT_EQ(sm, sn, sub_u, sub_ref);
}


MN_CASE( matview_whole_whole )
{
	test_mat_range<ref_matrix, whole, whole, M, N>();
}

MN_CASE( matview_whole_range )
{
	test_mat_range<ref_matrix, whole, range, M, N>();
}

MN_CASE( matview_whole_step )
{
	test_mat_range<ref_matrix, whole, step_range, M, N>();
}

MN_CASE( matview_range_whole )
{
	test_mat_range<ref_matrix, range, whole, M, N>();
}

MN_CASE( matview_range_range )
{
	test_mat_range<ref_matrix, range, range, M, N>();
}

MN_CASE( matview_range_step )
{
	test_mat_range<ref_matrix, range, step_range, M, N>();
}

MN_CASE( matview_step_whole )
{
	test_mat_range<ref_matrix, step_range, whole, M, N>();
}

MN_CASE( matview_step_range )
{
	test_mat_range<ref_matrix, step_range, range, M, N>();
}

MN_CASE( matview_step_step )
{
	test_mat_range<ref_matrix, step_range, step_range, M, N>();
}


MN_CASE( blocview_whole_whole )
{
	test_mat_range<ref_block, whole, whole, M, N>();
}

MN_CASE( blocview_whole_range )
{
	test_mat_range<ref_block, whole, range, M, N>();
}

MN_CASE( blocview_whole_step )
{
	test_mat_range<ref_block, whole, step_range, M, N>();
}

MN_CASE( blocview_range_whole )
{
	test_mat_range<ref_block, range, whole, M, N>();
}

MN_CASE( blocview_range_range )
{
	test_mat_range<ref_block, range, range, M, N>();
}

MN_CASE( blocview_range_step )
{
	test_mat_range<ref_block, range, step_range, M, N>();
}

MN_CASE( blocview_step_whole )
{
	test_mat_range<ref_block, step_range, whole, M, N>();
}

MN_CASE( blocview_step_range )
{
	test_mat_range<ref_block, step_range, range, M, N>();
}

MN_CASE( blocview_step_step )
{
	test_mat_range<ref_block, step_range, step_range, M, N>();
}


MN_CASE( gridview_whole_whole )
{
	test_mat_range<ref_grid, whole, whole, M, N>();
}

MN_CASE( gridview_whole_range )
{
	test_mat_range<ref_grid, whole, range, M, N>();
}

MN_CASE( gridview_whole_step )
{
	test_mat_range<ref_grid, whole, step_range, M, N>();
}

MN_CASE( gridview_range_whole )
{
	test_mat_range<ref_grid, range, whole, M, N>();
}

MN_CASE( gridview_range_range )
{
	test_mat_range<ref_grid, range, range, M, N>();
}

MN_CASE( gridview_range_step )
{
	test_mat_range<ref_grid, range, step_range, M, N>();
}

MN_CASE( gridview_step_whole )
{
	test_mat_range<ref_grid, step_range, whole, M, N>();
}

MN_CASE( gridview_step_range )
{
	test_mat_range<ref_grid, step_range, range, M, N>();
}

MN_CASE( gridview_step_step )
{
	test_mat_range<ref_grid, step_range, step_range, M, N>();
}



AUTO_TPACK( whole_whole_of_mat )
{
	ADD_MN_CASE_3X3( matview_whole_whole, DM, DN )
}

AUTO_TPACK( whole_range_of_mat )
{
	ADD_MN_CASE_3X3( matview_whole_range, DM, DN )
}

AUTO_TPACK( whole_step_of_mat )
{
	ADD_MN_CASE_3X3( matview_whole_step, DM, DN )
}

AUTO_TPACK( range_whole_of_mat )
{
	ADD_MN_CASE_3X3( matview_range_whole, DM, DN )
}

AUTO_TPACK( range_range_of_mat )
{
	ADD_MN_CASE_3X3( matview_range_range, DM, DN )
}

AUTO_TPACK( range_step_of_mat )
{
	ADD_MN_CASE_3X3( matview_range_step, DM, DN )
}

AUTO_TPACK( step_whole_of_mat )
{
	ADD_MN_CASE_3X3( matview_step_whole, DM, DN )
}

AUTO_TPACK( step_range_of_mat )
{
	ADD_MN_CASE_3X3( matview_step_range, DM, DN )
}

AUTO_TPACK( step_step_of_mat )
{
	ADD_MN_CASE_3X3( matview_step_step, DM, DN )
}


AUTO_TPACK( whole_whole_of_block )
{
	ADD_MN_CASE_3X3( blocview_whole_whole, DM, DN )
}

AUTO_TPACK( whole_range_of_block )
{
	ADD_MN_CASE_3X3( blocview_whole_range, DM, DN )
}

AUTO_TPACK( whole_step_of_block )
{
	ADD_MN_CASE_3X3( blocview_whole_step, DM, DN )
}

AUTO_TPACK( range_whole_of_block )
{
	ADD_MN_CASE_3X3( blocview_range_whole, DM, DN )
}

AUTO_TPACK( range_range_of_block )
{
	ADD_MN_CASE_3X3( blocview_range_range, DM, DN )
}

AUTO_TPACK( range_step_of_block )
{
	ADD_MN_CASE_3X3( blocview_range_step, DM, DN )
}

AUTO_TPACK( step_whole_of_block )
{
	ADD_MN_CASE_3X3( blocview_step_whole, DM, DN )
}

AUTO_TPACK( step_range_of_block )
{
	ADD_MN_CASE_3X3( blocview_step_range, DM, DN )
}

AUTO_TPACK( step_step_of_block )
{
	ADD_MN_CASE_3X3( blocview_step_step, DM, DN )
}


AUTO_TPACK( whole_whole_of_grid )
{
	ADD_MN_CASE_3X3( gridview_whole_whole, DM, DN )
}

AUTO_TPACK( whole_range_of_grid )
{
	ADD_MN_CASE_3X3( gridview_whole_range, DM, DN )
}

AUTO_TPACK( whole_step_of_grid )
{
	ADD_MN_CASE_3X3( gridview_whole_step, DM, DN )
}

AUTO_TPACK( range_whole_of_grid )
{
	ADD_MN_CASE_3X3( gridview_range_whole, DM, DN )
}

AUTO_TPACK( range_range_of_grid )
{
	ADD_MN_CASE_3X3( gridview_range_range, DM, DN )
}

AUTO_TPACK( range_step_of_grid )
{
	ADD_MN_CASE_3X3( gridview_range_step, DM, DN )
}

AUTO_TPACK( step_whole_of_grid )
{
	ADD_MN_CASE_3X3( gridview_step_whole, DM, DN )
}

AUTO_TPACK( step_range_of_grid )
{
	ADD_MN_CASE_3X3( gridview_step_range, DM, DN )
}

AUTO_TPACK( step_step_of_grid )
{
	ADD_MN_CASE_3X3( gridview_step_step, DM, DN )
}




