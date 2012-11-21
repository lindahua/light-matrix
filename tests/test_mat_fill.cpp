/**
 * @file test_mat_fill.cpp
 *
 * @brief Unit testing for matrix filling
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_block.h>
#include <light_mat/matrix/ref_grid.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;

const index_t cs_dst = 12;
const index_t rs_dst = 3;

const index_t LDim = 12;

// Auxiliary classes

template<template<typename T1, int R1, int C1> class SClassT, int M, int N> struct mat_maker;

template<int M, int N>
struct mat_maker<ref_matrix, M, N>
{
	static ref_matrix<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_matrix<double, M, N>(p, m, n);
	}
};

template<int M, int N>
struct mat_maker<ref_block, M, N>
{
	static ref_block<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_block<double, M, N>(p, m, n, cs_dst);
	}
};

template<int M, int N>
struct mat_maker<ref_grid, M, N>
{
	static ref_grid<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_grid<double, M, N>(p, m, n, rs_dst, cs_dst);
	}
};

template<class Mat, typename T>
bool verify_all_equal(const Mat& mat, const T& v)
{
	const index_t m = mat.nrows();
	const index_t n = mat.ncolumns();

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (mat(i, j) != v) return false;
		}
	}

	return true;
}


template<
	template<typename T2, int R2, int C2> class DClassT,
	int M, int N>
void test_matrix_zero()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> dst_mem(LDim * n, fill(-1.0));

	DClassT<double, M, N> dst = mat_maker<DClassT, M, N>::get_dst(dst_mem.ptr_data(), m, n);
	double v = 0;
	zero(dst);

	ASSERT_TRUE( verify_all_equal(dst, v) );
}


template<
	template<typename T2, int R2, int C2> class DClassT,
	int M, int N>
void test_matrix_fill()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> dst_mem(LDim * n, zero());

	DClassT<double, M, N> dst = mat_maker<DClassT, M, N>::get_dst(dst_mem.ptr_data(), m, n);
	double v = 12;
	fill(dst, v);

	ASSERT_TRUE( verify_all_equal(dst, v) );
}


MN_CASE( mat_zero, cont )
{
	test_matrix_zero<ref_matrix, M, N>();
}

MN_CASE( mat_zero, bloc )
{
	test_matrix_zero<ref_block, M, N>();
}

MN_CASE( mat_zero, grid )
{
	test_matrix_zero<ref_grid, M, N>();
}


MN_CASE( mat_fill, cont )
{
	test_matrix_fill<ref_matrix, M, N>();
}

MN_CASE( mat_fill, bloc )
{
	test_matrix_fill<ref_block, M, N>();
}

MN_CASE( mat_fill, grid )
{
	test_matrix_fill<ref_grid, M, N>();
}


BEGIN_TPACK( mat_zero_cont )
	ADD_MN_CASE_3X3( mat_zero, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_zero_bloc )
	ADD_MN_CASE_3X3( mat_zero, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_zero_grid )
	ADD_MN_CASE_3X3( mat_zero, grid, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_cont )
	ADD_MN_CASE_3X3( mat_fill, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_bloc )
	ADD_MN_CASE_3X3( mat_fill, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_fill_grid )
	ADD_MN_CASE_3X3( mat_fill, grid, 3, 4 )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_zero_cont )
	ADD_TPACK( mat_zero_bloc )
	ADD_TPACK( mat_zero_grid )

	ADD_TPACK( mat_fill_cont )
	ADD_TPACK( mat_fill_bloc )
	ADD_TPACK( mat_fill_grid )
END_MAIN_SUITE





