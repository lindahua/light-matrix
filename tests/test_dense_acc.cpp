/**
 * @file test_dense_acc.cpp
 *
 * Unit testing of accessors for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/dense_accessors.h>
#include <light_mat/matexpr/matrix_acc_eval.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary facilities

inline void fill_lin(dblock<double>& a)
{
	const index_t n = a.nelems();
	for (index_t i = 0; i < n; ++i)
	{
		a[i] = double(i+1);
	}
}

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



template<
	template<typename T1, int M1, int N1> class ClassT1,
	template<typename T2, int M2, int N2> class ClassT2,
	typename AccCate, typename KerCate, int M, int N>
void test_acc_eval()
{
	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	const index_t max_size1 = mat_maker<ClassT1, M, N>::max_size(m, n);
	const index_t max_size2 = mat_maker<ClassT2, M, N>::max_size(m, n);

	dblock<double> sblk(max_size1, zero());
	dblock<double> dblk(max_size2, zero());

	fill_lin(sblk);

	typedef typename mat_maker<ClassT1, M, N>::cmat_t smat_t;
	typedef typename mat_maker<ClassT2, M, N>::mat_t dmat_t;

	smat_t smat = mat_maker<ClassT1, M, N>::get_cmat(sblk.ptr_data(), m, n);
	dmat_t dmat = mat_maker<ClassT2, M, N>::get_mat(dblk.ptr_data(), m, n);

	typedef macc_setting<AccCate, KerCate, M, N> setting_t;

	macc_eval_scheme<smat_t, dmat_t, setting_t> scheme(smat, dmat);
	scheme.evaluate();

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


MN_CASE( linear_scalar, mat_to_mat )
{
	test_acc_eval<ref_matrix, ref_matrix,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, mat_to_bloc )
{
	test_acc_eval<ref_matrix, ref_block,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, mat_to_grid )
{
	test_acc_eval<ref_matrix, ref_grid,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, bloc_to_mat )
{
	test_acc_eval<ref_block, ref_matrix,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, bloc_to_bloc )
{
	test_acc_eval<ref_block, ref_block,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, bloc_to_grid )
{
	test_acc_eval<ref_block, ref_grid,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, grid_to_mat )
{
	test_acc_eval<ref_grid, ref_matrix,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, grid_to_bloc )
{
	test_acc_eval<ref_grid, ref_block,
		linear_macc, scalar_kernel_t, M, N>();
}

MN_CASE( linear_scalar, grid_to_grid )
{
	test_acc_eval<ref_grid, ref_grid,
		linear_macc, scalar_kernel_t, M, N>();
}


MN_CASE( percol_scalar, mat_to_mat )
{
	test_acc_eval<ref_matrix, ref_matrix,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, mat_to_bloc )
{
	test_acc_eval<ref_matrix, ref_block,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, mat_to_grid )
{
	test_acc_eval<ref_matrix, ref_grid,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, bloc_to_mat )
{
	test_acc_eval<ref_block, ref_matrix,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, bloc_to_bloc )
{
	test_acc_eval<ref_block, ref_block,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, bloc_to_grid )
{
	test_acc_eval<ref_block, ref_grid,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, grid_to_mat )
{
	test_acc_eval<ref_grid, ref_matrix,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, grid_to_bloc )
{
	test_acc_eval<ref_grid, ref_block,
		percol_macc, scalar_kernel_t, M, N>();
}

MN_CASE( percol_scalar, grid_to_grid )
{
	test_acc_eval<ref_grid, ref_grid,
		percol_macc, scalar_kernel_t, M, N>();
}



BEGIN_TPACK( linear_scalar_mat_to_mat )
	ADD_MN_CASE_3X3( linear_scalar, mat_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_mat_to_bloc )
	ADD_MN_CASE( linear_scalar, mat_to_bloc, 1, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_bloc, 0, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_bloc, DM, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_bloc, 1, 0 )
	ADD_MN_CASE( linear_scalar, mat_to_bloc, 1, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_mat_to_grid )
	ADD_MN_CASE( linear_scalar, mat_to_grid, 1, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_grid, 0, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_grid, DM, 1 )
	ADD_MN_CASE( linear_scalar, mat_to_grid, 1, 0 )
	ADD_MN_CASE( linear_scalar, mat_to_grid, 1, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_bloc_to_mat )
	ADD_MN_CASE_3X3( linear_scalar, bloc_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_bloc_to_bloc )
	ADD_MN_CASE( linear_scalar, bloc_to_bloc, 1, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_bloc, 0, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_bloc, DM, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_bloc, 1, 0 )
	ADD_MN_CASE( linear_scalar, bloc_to_bloc, 1, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_bloc_to_grid )
	ADD_MN_CASE( linear_scalar, bloc_to_grid, 1, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_grid, 0, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_grid, DM, 1 )
	ADD_MN_CASE( linear_scalar, bloc_to_grid, 1, 0 )
	ADD_MN_CASE( linear_scalar, bloc_to_grid, 1, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_grid_to_mat )
	ADD_MN_CASE_3X3( linear_scalar, grid_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_grid_to_bloc )
	ADD_MN_CASE( linear_scalar, grid_to_bloc, 1, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_bloc, 0, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_bloc, DM, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_bloc, 1, 0 )
	ADD_MN_CASE( linear_scalar, grid_to_bloc, 1, DN )
END_TPACK

BEGIN_TPACK( linear_scalar_grid_to_grid )
	ADD_MN_CASE( linear_scalar, grid_to_grid, 1, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_grid, 0, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_grid, DM, 1 )
	ADD_MN_CASE( linear_scalar, grid_to_grid, 1, 0 )
	ADD_MN_CASE( linear_scalar, grid_to_grid, 1, DN )
END_TPACK


BEGIN_TPACK( percol_scalar_mat_to_mat )
	ADD_MN_CASE_3X3( percol_scalar, mat_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_mat_to_bloc )
	ADD_MN_CASE_3X3( percol_scalar, mat_to_bloc, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_mat_to_grid )
	ADD_MN_CASE_3X3( percol_scalar, mat_to_grid, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_bloc_to_mat )
	ADD_MN_CASE_3X3( percol_scalar, bloc_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_bloc_to_bloc )
	ADD_MN_CASE_3X3( percol_scalar, bloc_to_bloc, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_bloc_to_grid )
	ADD_MN_CASE_3X3( percol_scalar, bloc_to_grid, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_grid_to_mat )
	ADD_MN_CASE_3X3( percol_scalar, grid_to_mat, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_grid_to_bloc )
	ADD_MN_CASE_3X3( percol_scalar, grid_to_bloc, DM, DN )
END_TPACK

BEGIN_TPACK( percol_scalar_grid_to_grid )
	ADD_MN_CASE_3X3( percol_scalar, grid_to_grid, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( linear_scalar_mat_to_mat )
	ADD_TPACK( linear_scalar_mat_to_bloc )
	ADD_TPACK( linear_scalar_mat_to_grid )
	ADD_TPACK( linear_scalar_bloc_to_mat )
	ADD_TPACK( linear_scalar_bloc_to_bloc )
	ADD_TPACK( linear_scalar_bloc_to_grid )
	ADD_TPACK( linear_scalar_grid_to_mat )
	ADD_TPACK( linear_scalar_grid_to_bloc )
	ADD_TPACK( linear_scalar_grid_to_grid )

	ADD_TPACK( percol_scalar_mat_to_mat )
	ADD_TPACK( percol_scalar_mat_to_bloc )
	ADD_TPACK( percol_scalar_mat_to_grid )
	ADD_TPACK( percol_scalar_bloc_to_mat )
	ADD_TPACK( percol_scalar_bloc_to_bloc )
	ADD_TPACK( percol_scalar_bloc_to_grid )
	ADD_TPACK( percol_scalar_grid_to_mat )
	ADD_TPACK( percol_scalar_grid_to_bloc )
	ADD_TPACK( percol_scalar_grid_to_grid )
END_MAIN_SUITE




