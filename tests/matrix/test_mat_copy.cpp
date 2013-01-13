/**
 * @file test_mat_copy.cpp
 *
 * Unit testing of matrix copy
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_block.h>
#include <light_mat/matrix/ref_grid.h>
#include <light_mat/common/block.h>

using namespace lmat;
using namespace lmat::test;


// Auxiliary classes


template<typename STag, typename DTag, int M, int N>
void test_matrix_copy()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<STag, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	mat_host<STag, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	copy(smat, dmat);

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


template<typename DTag, int M, int N>
void test_matrix_import()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<cont, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<DTag, double, M, N>::mat_t dmat_t;

	mat_host<cont, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<DTag, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	copy(smat.ptr_data(), dmat);

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


template<typename STag, int M, int N>
void test_matrix_export()
{
	const index_t m = M == 0 ? 3 : M;
	const index_t n = N == 0 ? 4 : N;

	typedef typename mat_host<STag, double, M, N>::cmat_t smat_t;
	typedef typename mat_host<cont, double, M, N>::mat_t dmat_t;

	mat_host<STag, double, M, N> src(m, n);
	src.fill_lin();
	mat_host<cont, double, M, N> dst(m, n);

	smat_t smat = src.get_cmat();
	dmat_t dmat = dst.get_mat();

	copy(smat, dmat.ptr_data());

	ASSERT_MAT_EQ(m, n, smat, dmat);
}


MN_CASE( mat_copy_cont_to_cont )
{
	test_matrix_copy<cont, cont, M, N>();
}

MN_CASE( mat_copy_cont_to_bloc )
{
	test_matrix_copy<cont, bloc, M, N>();
}

MN_CASE( mat_copy_cont_to_grid )
{
	test_matrix_copy<cont, grid, M, N>();
}

MN_CASE( mat_copy_bloc_to_cont )
{
	test_matrix_copy<bloc, cont, M, N>();
}

MN_CASE( mat_copy_bloc_to_bloc )
{
	test_matrix_copy<bloc, bloc, M, N>();
}

MN_CASE( mat_copy_bloc_to_grid )
{
	test_matrix_copy<bloc, grid, M, N>();
}

MN_CASE( mat_copy_grid_to_cont )
{
	test_matrix_copy<grid, cont, M, N>();
}

MN_CASE( mat_copy_grid_to_bloc )
{
	test_matrix_copy<grid, bloc, M, N>();
}

MN_CASE( mat_copy_grid_to_grid )
{
	test_matrix_copy<grid, grid, M, N>();
}


MN_CASE( mat_import_cont )
{
	test_matrix_import<cont, M, N>();
}

MN_CASE( mat_import_bloc )
{
	test_matrix_import<bloc, M, N>();
}

MN_CASE( mat_import_grid )
{
	test_matrix_import<grid, M, N>();
}


MN_CASE( mat_export_cont )
{
	test_matrix_export<cont, M, N>();
}

MN_CASE( mat_export_bloc )
{
	test_matrix_export<bloc, M, N>();
}

MN_CASE( mat_export_grid )
{
	test_matrix_export<grid, M, N>();
}


template<class S, class D>
void safe_copy_triu(index_t m, index_t n, const S& smat, D& dmat, index_t k)
{
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (j - i >= k) dmat(i, j) = smat(i, j);
		}
	}
}

template<class S, class D>
void safe_copy_tril(index_t m, index_t n, const S& smat, D& dmat, index_t k)
{
	for (index_t j = 0; j < n; ++j)
	{
		for (index_t i = 0; i < m; ++i)
		{
			if (j - i <= k) dmat(i, j) = smat(i, j);
		}
	}
}



SIMPLE_CASE( mat_copy_triu )
{
	index_t LDimA = 9;
	index_t LDimB = 10;

	index_t max_m = 8;
	index_t max_n = 8;

	dense_matrix<double> A(LDimA, max_n);
	for (index_t i = 0; i < A.nelems(); ++i) A[i] = double(i+1);

	dense_matrix<double> B(LDimB, max_n);

	for (index_t m = 1; m < max_m; ++m)
	{
		for (index_t n = 1; n < max_n; ++n)
		{
			cref_block<double> a(A.ptr_data(), m, n, LDimA);
			ref_block<double> b(B.ptr_data(), m, n, LDimB);

			dense_matrix<double> r(m, n);

			index_t max_k = m < n ? m - 1 : n - 1;
			for (index_t k = -max_k; k <= max_k; ++k)
			{
				zero_vec(B.nelems(), B.ptr_data());
				zero_vec(r.nelems(), r.ptr_data());

				copy_triu(a, b, k);
				safe_copy_triu(m, n, a, r, k);

				ASSERT_MAT_EQ(m, n, b, r);
			}
		}
	}
}


SIMPLE_CASE( mat_copy_tril )
{
	index_t LDimA = 9;
	index_t LDimB = 10;

	index_t max_m = 8;
	index_t max_n = 8;

	dense_matrix<double> A(LDimA, max_n);
	for (index_t i = 0; i < A.nelems(); ++i) A[i] = double(i+1);

	dense_matrix<double> B(LDimB, max_n);

	for (index_t m = 1; m < max_m; ++m)
	{
		for (index_t n = 1; n < max_n; ++n)
		{
			cref_block<double> a(A.ptr_data(), m, n, LDimA);
			ref_block<double> b(B.ptr_data(), m, n, LDimB);

			dense_matrix<double> r(m, n);

			index_t max_k = m < n ? m - 1 : n - 1;
			for (index_t k = -max_k; k <= max_k; ++k)
			{
				zero_vec(B.nelems(), B.ptr_data());
				zero_vec(r.nelems(), r.ptr_data());

				copy_tril(a, b, k);
				safe_copy_tril(m, n, a, r, k);

				ASSERT_MAT_EQ(m, n, b, r);
			}
		}
	}
}


LTEST_INIT_AUTOSUITE

AUTO_TPACK( mat_copy )
{
	ADD_MN_CASE_3X3( mat_copy_cont_to_cont, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_cont_to_bloc, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_cont_to_grid, 3, 4 )

	ADD_MN_CASE_3X3( mat_copy_bloc_to_cont, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_bloc_to_bloc, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_bloc_to_grid, 3, 4 )

	ADD_MN_CASE_3X3( mat_copy_grid_to_cont, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_grid_to_bloc, 3, 4 )
	ADD_MN_CASE_3X3( mat_copy_grid_to_grid, 3, 4 )
}

AUTO_TPACK( mat_import )
{
	ADD_MN_CASE_3X3( mat_import_cont, 3, 4 )
	ADD_MN_CASE_3X3( mat_import_bloc, 3, 4 )
	ADD_MN_CASE_3X3( mat_import_grid, 3, 4 )
}

AUTO_TPACK( mat_export )
{
	ADD_MN_CASE_3X3( mat_export_cont, 3, 4 )
	ADD_MN_CASE_3X3( mat_export_bloc, 3, 4 )
	ADD_MN_CASE_3X3( mat_export_grid, 3, 4 )
}

AUTO_TPACK( mat_copy_tri )
{
	ADD_SIMPLE_CASE( mat_copy_triu )
	ADD_SIMPLE_CASE( mat_copy_tril )
}

