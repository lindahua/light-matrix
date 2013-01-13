/**
 * @file test_mat_compare.cpp
 *
 * Unit testing for matrix comparison
 *
 * @author Dahua Lin
 */


#include "../test_base.h"
#include <light_mat/mateval/mat_compare.h>

using namespace lmat;
using namespace lmat::test;

using namespace lmat;
using namespace lmat::test;

const index_t cs_a = 11;
const index_t cs_b = 12;
const index_t rs_a = 2;
const index_t rs_b = 3;

const index_t LDim = 12;

// Auxiliary classes

template<template<typename T1, int R1, int C1> class SClassT, int M, int N> struct mat_maker;

template<int M, int N>
struct mat_maker<ref_matrix, M, N>
{
	static ref_matrix<double, M, N> get_a(double *p, index_t m, index_t n)
	{
		return ref_matrix<double, M, N>(p, m, n);
	}

	static ref_matrix<double, M, N> get_b(double *p, index_t m, index_t n)
	{
		return ref_matrix<double, M, N>(p, m, n);
	}
};

template<int M, int N>
struct mat_maker<ref_block, M, N>
{
	static ref_block<double, M, N> get_a(double *p, index_t m, index_t n)
	{
		return ref_block<double, M, N>(p, m, n, cs_a);
	}

	static ref_block<double, M, N> get_b(double *p, index_t m, index_t n)
	{
		return ref_block<double, M, N>(p, m, n, cs_b);
	}
};

template<int M, int N>
struct mat_maker<ref_grid, M, N>
{
	static ref_grid<double, M, N> get_a(double *p, index_t m, index_t n)
	{
		return ref_grid<double, M, N>(p, m, n, rs_a, cs_a);
	}

	static ref_grid<double, M, N> get_b(double *p, index_t m, index_t n)
	{
		return ref_grid<double, M, N>(p, m, n, rs_b, cs_b);
	}
};



template<
	template<typename T1, int R1, int C1> class AClassT,
	template<typename T2, int R2, int C2> class BClassT,
	int M, int N>
void test_matrix_equal()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> a_mem(LDim * n, zero());
	dblock<double> b_mem(LDim * n, zero());

	AClassT<double, M, N> a = mat_maker<AClassT, M, N>::get_a(a_mem.ptr_data(), m, n);
	BClassT<double, M, N> b = mat_maker<BClassT, M, N>::get_b(b_mem.ptr_data(), m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			a(i, j) = b(i, j) = double(2 * i + 3 * j + 1);
		}
	}

	ASSERT_TRUE( is_equal(a, b) );

	b(m-1, n-1) += 1.0;

	ASSERT_FALSE( is_equal(a, b) );
}


template<
	template<typename T1, int R1, int C1> class AClassT,
	template<typename T2, int R2, int C2> class BClassT,
	int M, int N>
void test_matrix_approx()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	dblock<double> a_mem(LDim * n, zero());
	dblock<double> b_mem(LDim * n, zero());

	AClassT<double, M, N> a = mat_maker<AClassT, M, N>::get_a(a_mem.ptr_data(), m, n);
	BClassT<double, M, N> b = mat_maker<BClassT, M, N>::get_b(b_mem.ptr_data(), m, n);

	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			a(i, j) = b(i, j) = double(2 * i + 3 * j + 1);
		}
	}

	ASSERT_TRUE( is_approx(a, b, 1.0e-3) );

	b(m-1, n-1) += 1.0;

	ASSERT_FALSE( is_approx(a, b, 1.0e-3) );
	ASSERT_TRUE( is_approx(a, b, 1.2) );
}


MN_CASE( mat_equal_cont_to_cont )
{
	test_matrix_equal<ref_matrix, ref_matrix, M, N>();
}

MN_CASE( mat_equal_cont_to_bloc )
{
	test_matrix_equal<ref_matrix, ref_block, M, N>();
}

MN_CASE( mat_equal_cont_to_grid )
{
	test_matrix_equal<ref_matrix, ref_grid, M, N>();
}

MN_CASE( mat_equal_bloc_to_cont )
{
	test_matrix_equal<ref_block, ref_matrix, M, N>();
}

MN_CASE( mat_equal_bloc_to_bloc )
{
	test_matrix_equal<ref_block, ref_block, M, N>();
}

MN_CASE( mat_equal_bloc_to_grid )
{
	test_matrix_equal<ref_block, ref_grid, M, N>();
}

MN_CASE( mat_equal_grid_to_cont )
{
	test_matrix_equal<ref_grid, ref_matrix, M, N>();
}

MN_CASE( mat_equal_grid_to_bloc )
{
	test_matrix_equal<ref_grid, ref_block, M, N>();
}

MN_CASE( mat_equal_grid_to_grid )
{
	test_matrix_equal<ref_grid, ref_grid, M, N>();
}


MN_CASE( mat_approx_cont_to_cont )
{
	test_matrix_approx<ref_matrix, ref_matrix, M, N>();
}

MN_CASE( mat_approx_cont_to_bloc )
{
	test_matrix_approx<ref_matrix, ref_block, M, N>();
}

MN_CASE( mat_approx_cont_to_grid )
{
	test_matrix_approx<ref_matrix, ref_grid, M, N>();
}

MN_CASE( mat_approx_bloc_to_cont )
{
	test_matrix_approx<ref_block, ref_matrix, M, N>();
}

MN_CASE( mat_approx_bloc_to_bloc )
{
	test_matrix_approx<ref_block, ref_block, M, N>();
}

MN_CASE( mat_approx_bloc_to_grid )
{
	test_matrix_approx<ref_block, ref_grid, M, N>();
}

MN_CASE( mat_approx_grid_to_cont )
{
	test_matrix_approx<ref_grid, ref_matrix, M, N>();
}

MN_CASE( mat_approx_grid_to_bloc )
{
	test_matrix_approx<ref_grid, ref_block, M, N>();
}

MN_CASE( mat_approx_grid_to_grid )
{
	test_matrix_approx<ref_grid, ref_grid, M, N>();
}



AUTO_TPACK( mat_equal_cc )
{
	ADD_MN_CASE_3X3( mat_equal_cont_to_cont, 3, 4 );
}

AUTO_TPACK( mat_equal_cb )
{
	ADD_MN_CASE_3X3( mat_equal_cont_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_equal_cg )
{
	ADD_MN_CASE_3X3( mat_equal_cont_to_grid, 3, 4 );
}

AUTO_TPACK( mat_equal_bc )
{
	ADD_MN_CASE_3X3( mat_equal_bloc_to_cont, 3, 4 );
}

AUTO_TPACK( mat_equal_bb )
{
	ADD_MN_CASE_3X3( mat_equal_bloc_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_equal_bg )
{
	ADD_MN_CASE_3X3( mat_equal_bloc_to_grid, 3, 4 );
}

AUTO_TPACK( mat_equal_gc )
{
	ADD_MN_CASE_3X3( mat_equal_grid_to_cont, 3, 4 );
}

AUTO_TPACK( mat_equal_gb )
{
	ADD_MN_CASE_3X3( mat_equal_grid_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_equal_gg )
{
	ADD_MN_CASE_3X3( mat_equal_grid_to_grid, 3, 4 );
}


AUTO_TPACK( mat_approx_cc )
{
	ADD_MN_CASE_3X3( mat_approx_cont_to_cont, 3, 4 );
}

AUTO_TPACK( mat_approx_cb )
{
	ADD_MN_CASE_3X3( mat_approx_cont_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_approx_cg )
{
	ADD_MN_CASE_3X3( mat_approx_cont_to_grid, 3, 4 );
}

AUTO_TPACK( mat_approx_bc )
{
	ADD_MN_CASE_3X3( mat_approx_bloc_to_cont, 3, 4 );
}

AUTO_TPACK( mat_approx_bb )
{
	ADD_MN_CASE_3X3( mat_approx_bloc_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_approx_bg )
{
	ADD_MN_CASE_3X3( mat_approx_bloc_to_grid, 3, 4 );
}

AUTO_TPACK( mat_approx_gc )
{
	ADD_MN_CASE_3X3( mat_approx_grid_to_cont, 3, 4 );
}

AUTO_TPACK( mat_approx_gb )
{
	ADD_MN_CASE_3X3( mat_approx_grid_to_bloc, 3, 4 );
}

AUTO_TPACK( mat_approx_gg )
{
	ADD_MN_CASE_3X3( mat_approx_grid_to_grid, 3, 4 );
}


