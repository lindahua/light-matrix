/**
 * @file test_mat_copy.cpp
 *
 * Unit testing of matrix copy
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

const index_t cs_src = 11;
const index_t cs_dst = 12;
const index_t rs_src = 2;
const index_t rs_dst = 3;

const index_t LDim = 12;

// Auxiliary classes

template<template<typename T1, int R1, int C1> class SClassT, int M, int N> struct mat_maker;

template<int M, int N>
struct mat_maker<ref_matrix, M, N>
{
	static ref_matrix<double, M, N> get_src(double *p, index_t m, index_t n)
	{
		return ref_matrix<double, M, N>(p, m, n);
	}

	static ref_matrix<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_matrix<double, M, N>(p, m, n);
	}
};

template<int M, int N>
struct mat_maker<ref_block, M, N>
{
	static ref_block<double, M, N> get_src(double *p, index_t m, index_t n)
	{
		return ref_block<double, M, N>(p, m, n, cs_src);
	}

	static ref_block<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_block<double, M, N>(p, m, n, cs_dst);
	}
};

template<int M, int N>
struct mat_maker<ref_grid, M, N>
{
	static ref_grid<double, M, N> get_src(double *p, index_t m, index_t n)
	{
		return ref_grid<double, M, N>(p, m, n, rs_src, cs_src);
	}

	static ref_grid<double, M, N> get_dst(double *p, index_t m, index_t n)
	{
		return ref_grid<double, M, N>(p, m, n, rs_dst, cs_dst);
	}
};



template<
	template<typename T1, int R1, int C1> class SClassT,
	template<typename T2, int R2, int C2> class DClassT,
	int M, int N>
void test_matrix_copy()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> src_mem(LDim * n, zero());
	dblock<double> dst_mem(LDim * n, zero());

	SClassT<double, M, N> src = mat_maker<SClassT, M, N>::get_src(src_mem.ptr_data(), m, n);
	DClassT<double, M, N> dst = mat_maker<DClassT, M, N>::get_dst(dst_mem.ptr_data(), m, n);

	for (index_t i = 0; i < LDim * n; ++i) src_mem[i] = double(2 * i + 1);
	copy(src, dst);

	ASSERT_MAT_EQ(m, n, src, dst);
}


template<
	template<typename T2, int R2, int C2> class DClassT,
	int M, int N>
void test_matrix_import()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> src_mem(LDim * n, zero());
	dblock<double> dst_mem(LDim * n, zero());

	DClassT<double, M, N> dst = mat_maker<DClassT, M, N>::get_dst(dst_mem.ptr_data(), m, n);

	for (index_t i = 0; i < LDim * n; ++i) src_mem[i] = double(2 * i + 1);
	copy(src_mem.ptr_data(), dst);

	ref_matrix<double, M, N> src(src_mem.ptr_data(), m, n);
	ASSERT_MAT_EQ(m, n, src, dst);
}


template<
	template<typename T1, int R1, int C1> class SClassT,
	int M, int N>
void test_matrix_export()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> src_mem(LDim * n, zero());
	dblock<double> dst_mem(LDim * n, zero());

	SClassT<double, M, N> src = mat_maker<SClassT, M, N>::get_src(src_mem.ptr_data(), m, n);

	for (index_t i = 0; i < LDim * n; ++i) src_mem[i] = double(2 * i + 1);
	copy(src, dst_mem.ptr_data());

	ref_matrix<double, M, N> dst(dst_mem.ptr_data(), m, n);
	ASSERT_MAT_EQ(m, n, src, dst);
}



MN_CASE( mat_copy, cont_to_cont )
{
	test_matrix_copy<ref_matrix, ref_matrix, M, N>();
}

MN_CASE( mat_copy, cont_to_bloc )
{
	test_matrix_copy<ref_matrix, ref_block, M, N>();
}

MN_CASE( mat_copy, cont_to_grid )
{
	test_matrix_copy<ref_matrix, ref_grid, M, N>();
}

MN_CASE( mat_copy, bloc_to_cont )
{
	test_matrix_copy<ref_block, ref_matrix, M, N>();
}

MN_CASE( mat_copy, bloc_to_bloc )
{
	test_matrix_copy<ref_block, ref_block, M, N>();
}

MN_CASE( mat_copy, bloc_to_grid )
{
	test_matrix_copy<ref_block, ref_grid, M, N>();
}

MN_CASE( mat_copy, grid_to_cont )
{
	test_matrix_copy<ref_grid, ref_matrix, M, N>();
}

MN_CASE( mat_copy, grid_to_bloc )
{
	test_matrix_copy<ref_grid, ref_block, M, N>();
}

MN_CASE( mat_copy, grid_to_grid )
{
	test_matrix_copy<ref_grid, ref_grid, M, N>();
}


MN_CASE( mat_import, cont )
{
	test_matrix_import<ref_matrix, M, N>();
}

MN_CASE( mat_import, bloc )
{
	test_matrix_import<ref_block, M, N>();
}

MN_CASE( mat_import, grid )
{
	test_matrix_import<ref_grid, M, N>();
}


MN_CASE( mat_export, cont )
{
	test_matrix_export<ref_matrix, M, N>();
}

MN_CASE( mat_export, bloc )
{
	test_matrix_export<ref_block, M, N>();
}

MN_CASE( mat_export, grid )
{
	test_matrix_export<ref_grid, M, N>();
}



BEGIN_TPACK( mat_copy_cc )
	ADD_MN_CASE_3X3( mat_copy, cont_to_cont, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_cb )
	ADD_MN_CASE_3X3( mat_copy, cont_to_bloc, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_cg )
	ADD_MN_CASE_3X3( mat_copy, cont_to_grid, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_bc )
	ADD_MN_CASE_3X3( mat_copy, bloc_to_cont, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_bb )
	ADD_MN_CASE_3X3( mat_copy, bloc_to_bloc, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_bg )
	ADD_MN_CASE_3X3( mat_copy, bloc_to_grid, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_gc )
	ADD_MN_CASE_3X3( mat_copy, grid_to_cont, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_gb )
	ADD_MN_CASE_3X3( mat_copy, grid_to_bloc, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_copy_gg )
	ADD_MN_CASE_3X3( mat_copy, grid_to_grid, 3, 4 );
END_TPACK

BEGIN_TPACK( mat_import_cont )
	ADD_MN_CASE_3X3( mat_import, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_import_bloc )
	ADD_MN_CASE_3X3( mat_import, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_import_grid )
	ADD_MN_CASE_3X3( mat_import, grid, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_export_cont )
	ADD_MN_CASE_3X3( mat_export, cont, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_export_bloc )
	ADD_MN_CASE_3X3( mat_export, bloc, 3, 4 )
END_TPACK

BEGIN_TPACK( mat_export_grid )
	ADD_MN_CASE_3X3( mat_export, grid, 3, 4 )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( mat_copy_cc )
	ADD_TPACK( mat_copy_cb )
	ADD_TPACK( mat_copy_cg )
	ADD_TPACK( mat_copy_bc )
	ADD_TPACK( mat_copy_bb )
	ADD_TPACK( mat_copy_bg )
	ADD_TPACK( mat_copy_gc )
	ADD_TPACK( mat_copy_gb )
	ADD_TPACK( mat_copy_gg )

	ADD_TPACK( mat_import_cont )
	ADD_TPACK( mat_import_bloc )
	ADD_TPACK( mat_import_grid )

	ADD_TPACK( mat_export_cont )
	ADD_TPACK( mat_export_bloc )
	ADD_TPACK( mat_export_grid )
END_MAIN_SUITE



