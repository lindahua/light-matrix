/**
 * @file test_mat_fold.cpp
 *
 * @brief Test matrix folding
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#define DEFAULT_M_VALUE 9
#define DEFAULT_N_VALUE 8

#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_fold.h>
#include <limits>

using namespace lmat;
using namespace lmat::test;


struct sum_tt
{
	typedef sum_kernel<double> kernel_type;

	static double tol()
	{
		return 1.0e-12;
	}

	template<class A>
	static double eval(const IRegularMatrix<A, double>& a)
	{
		double r(0);
		for (index_t j = 0; j < a.ncolumns(); ++j)
		{
			for (index_t i = 0; i < a.nrows(); ++i)
				r += a(i, j);
		}

		return r;
	}
};

struct max_tt
{
	typedef maximum_kernel<double> kernel_type;

	static double tol()
	{
		return 1.0e-16;
	}

	template<class A>
	static double eval(const IRegularMatrix<A, double>& a)
	{
		double r(-std::numeric_limits<double>::infinity());
		for (index_t j = 0; j < a.ncolumns(); ++j)
		{
			for (index_t i = 0; i < a.nrows(); ++i)
				if (a(i, j) > r) r = a(i, j);
		}

		return r;
	}
};

struct min_tt
{
	typedef minimum_kernel<double> kernel_type;

	static double tol()
	{
		return 1.0e-16;
	}

	template<class A>
	static double eval(const IRegularMatrix<A, double>& a)
	{
		double r(std::numeric_limits<double>::infinity());
		for (index_t j = 0; j < a.ncolumns(); ++j)
		{
			for (index_t i = 0; i < a.nrows(); ++i)
				if (a(i, j) < r) r = a(i, j);
		}

		return r;
	}
};


template<index_t CM, index_t CN>
inline bool my_use_linear(cont, matrix_shape<CM, CN>)
{
	return true;
}

template<index_t CM, index_t CN>
inline bool my_use_linear(bloc, matrix_shape<CM, CN>)
{
	return CN == 1 || CM == 1;
}

template<index_t CM, index_t CN>
inline bool my_use_linear(grid, matrix_shape<CM, CN>)
{
	return CN == 1 || CM == 1;
}

template<index_t CM, index_t CN>
inline bool my_use_simd(cont, matrix_shape<CM, CN>)
{
	return (CM * CN) % 4 == 0;
}

template<index_t CM, index_t CN>
inline bool my_use_simd(bloc, matrix_shape<CM, CN>)
{
	return CM % 4 == 0;
}

template<index_t CM, index_t CN>
inline bool my_use_simd(grid, matrix_shape<CM, CN>)
{
	return false;
}



template<class KTT, typename Acc, typename U, class MTag, index_t CM, index_t CN>
void test_folder_x()
{
	typedef typename KTT::kernel_type kernel_t;
	kernel_t fker;
	typedef typename kernel_t::value_type VT;
	typedef typename mat_host<MTag, VT, CM, CN>::cmat_t smat_t;

	const index_t m = CM == 0 ? DM : CM;
	const index_t n = CN == 0 ? DN : CN;

	mat_host<MTag, double, CM, CN> s(m, n);
	s.fill_rand();
	smat_t smat = s.get_cmat();

	VT r0 = KTT::eval(smat);

	VT tol = KTT::tol();

	VT r1 = fold(fker).eval(macc_<Acc, U>(), m, n, in_(smat));
	ASSERT_APPROX(r1, r0, tol);

	matrix_shape<CM, CN> shape(m, n);
	VT r2 = fold(fker).eval(macc_<Acc, U>(), shape, in_(smat));
	ASSERT_APPROX(r2, r0, tol);
}


template<class KTT, class MTag, index_t CM, index_t CN>
void test_folder()
{
	typedef typename KTT::kernel_type kernel_t;
	kernel_t fker;
	typedef typename kernel_t::value_type VT;
	typedef typename mat_host<MTag, VT, CM, CN>::cmat_t smat_t;

	const index_t m = CM == 0 ? DM : CM;
	const index_t n = CN == 0 ? DN : CN;

	mat_host<MTag, double, CM, CN> s(m, n);
	s.fill_rand();
	smat_t smat = s.get_cmat();

	matrix_shape<CM, CN> shape(m, n);

	bool expect_use_linear = my_use_linear(MTag(), shape);
	bool expect_use_simd = my_use_simd(MTag(), shape);

	typename internal::fold_policy<kernel_t, matrix_shape<CM, CN>, smat_t>::type policy;
	ASSERT_EQ( use_linear_acc(policy), expect_use_linear );
	ASSERT_EQ( use_simd(policy), expect_use_simd );

	VT r0 = KTT::eval(smat);

	VT tol = KTT::tol();
	VT r1 = fold(fker)(m, n, in_(smat));
	ASSERT_APPROX(r1, r0, tol);

	VT r2 = fold(fker)(shape, in_(smat));
	ASSERT_APPROX(r2, r0, tol);
}




// specific cases

// sum

MN_CASE( sum_linear_scalar_cont )
{
	test_folder_x<sum_tt, linear_, scalar_, cont, M, N>();
}

MN_CASE( sum_linear_scalar_bloc )
{
	test_folder_x<sum_tt, linear_, scalar_, bloc, M, N>();
}

MN_CASE( sum_linear_scalar_grid )
{
	test_folder_x<sum_tt, linear_, scalar_, grid, M, N>();
}

MN_CASE( sum_linear_sse_cont )
{
	test_folder_x<sum_tt, linear_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( sum_linear_sse_bloc )
{
	test_folder_x<sum_tt, linear_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( sum_linear_sse_grid )
{
	test_folder_x<sum_tt, linear_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( sum_linear_avx_cont )
{
	test_folder_x<sum_tt, linear_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( sum_linear_avx_bloc )
{
	test_folder_x<sum_tt, linear_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( sum_linear_avx_grid )
{
	test_folder_x<sum_tt, linear_, simd_<avx_t>, grid, M, N>();
}

#endif


MN_CASE( sum_percol_scalar_cont )
{
	test_folder_x<sum_tt, percol_, scalar_, cont, M, N>();
}

MN_CASE( sum_percol_scalar_bloc )
{
	test_folder_x<sum_tt, percol_, scalar_, bloc, M, N>();
}

MN_CASE( sum_percol_scalar_grid )
{
	test_folder_x<sum_tt, percol_, scalar_, grid, M, N>();
}

MN_CASE( sum_percol_sse_cont )
{
	test_folder_x<sum_tt, percol_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( sum_percol_sse_bloc )
{
	test_folder_x<sum_tt, percol_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( sum_percol_sse_grid )
{
	test_folder_x<sum_tt, percol_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( sum_percol_avx_cont )
{
	test_folder_x<sum_tt, percol_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( sum_percol_avx_bloc )
{
	test_folder_x<sum_tt, percol_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( sum_percol_avx_grid )
{
	test_folder_x<sum_tt, percol_, simd_<avx_t>, grid, M, N>();
}

#endif

MN_CASE( sum_auto_cont )
{
	test_folder<sum_tt, cont, M, N>();
}

MN_CASE( sum_auto_bloc )
{
	test_folder<sum_tt, bloc, M, N>();
}

MN_CASE( sum_auto_grid )
{
	test_folder<sum_tt, grid, M, N>();
}


// max

MN_CASE( max_linear_scalar_cont )
{
	test_folder_x<max_tt, linear_, scalar_, cont, M, N>();
}

MN_CASE( max_linear_scalar_bloc )
{
	test_folder_x<max_tt, linear_, scalar_, bloc, M, N>();
}

MN_CASE( max_linear_scalar_grid )
{
	test_folder_x<max_tt, linear_, scalar_, grid, M, N>();
}

MN_CASE( max_linear_sse_cont )
{
	test_folder_x<max_tt, linear_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( max_linear_sse_bloc )
{
	test_folder_x<max_tt, linear_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( max_linear_sse_grid )
{
	test_folder_x<max_tt, linear_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( max_linear_avx_cont )
{
	test_folder_x<max_tt, linear_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( max_linear_avx_bloc )
{
	test_folder_x<max_tt, linear_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( max_linear_avx_grid )
{
	test_folder_x<max_tt, linear_, simd_<avx_t>, grid, M, N>();
}

#endif


MN_CASE( max_percol_scalar_cont )
{
	test_folder_x<max_tt, percol_, scalar_, cont, M, N>();
}

MN_CASE( max_percol_scalar_bloc )
{
	test_folder_x<max_tt, percol_, scalar_, bloc, M, N>();
}

MN_CASE( max_percol_scalar_grid )
{
	test_folder_x<max_tt, percol_, scalar_, grid, M, N>();
}

MN_CASE( max_percol_sse_cont )
{
	test_folder_x<max_tt, percol_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( max_percol_sse_bloc )
{
	test_folder_x<max_tt, percol_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( max_percol_sse_grid )
{
	test_folder_x<max_tt, percol_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( max_percol_avx_cont )
{
	test_folder_x<max_tt, percol_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( max_percol_avx_bloc )
{
	test_folder_x<max_tt, percol_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( max_percol_avx_grid )
{
	test_folder_x<max_tt, percol_, simd_<avx_t>, grid, M, N>();
}

#endif

MN_CASE( max_auto_cont )
{
	test_folder<max_tt, cont, M, N>();
}

MN_CASE( max_auto_bloc )
{
	test_folder<max_tt, bloc, M, N>();
}

MN_CASE( max_auto_grid )
{
	test_folder<max_tt, grid, M, N>();
}


// min

MN_CASE( min_linear_scalar_cont )
{
	test_folder_x<min_tt, linear_, scalar_, cont, M, N>();
}

MN_CASE( min_linear_scalar_bloc )
{
	test_folder_x<min_tt, linear_, scalar_, bloc, M, N>();
}

MN_CASE( min_linear_scalar_grid )
{
	test_folder_x<min_tt, linear_, scalar_, grid, M, N>();
}

MN_CASE( min_linear_sse_cont )
{
	test_folder_x<min_tt, linear_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( min_linear_sse_bloc )
{
	test_folder_x<min_tt, linear_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( min_linear_sse_grid )
{
	test_folder_x<min_tt, linear_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( min_linear_avx_cont )
{
	test_folder_x<min_tt, linear_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( min_linear_avx_bloc )
{
	test_folder_x<min_tt, linear_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( min_linear_avx_grid )
{
	test_folder_x<min_tt, linear_, simd_<avx_t>, grid, M, N>();
}

#endif


MN_CASE( min_percol_scalar_cont )
{
	test_folder_x<min_tt, percol_, scalar_, cont, M, N>();
}

MN_CASE( min_percol_scalar_bloc )
{
	test_folder_x<min_tt, percol_, scalar_, bloc, M, N>();
}

MN_CASE( min_percol_scalar_grid )
{
	test_folder_x<min_tt, percol_, scalar_, grid, M, N>();
}

MN_CASE( min_percol_sse_cont )
{
	test_folder_x<min_tt, percol_, simd_<sse_t>, cont, M, N>();
}

MN_CASE( min_percol_sse_bloc )
{
	test_folder_x<min_tt, percol_, simd_<sse_t>, bloc, M, N>();
}

MN_CASE( min_percol_sse_grid )
{
	test_folder_x<min_tt, percol_, simd_<sse_t>, grid, M, N>();
}

#ifdef LMAT_HAS_AVX

MN_CASE( min_percol_avx_cont )
{
	test_folder_x<min_tt, percol_, simd_<avx_t>, cont, M, N>();
}

MN_CASE( min_percol_avx_bloc )
{
	test_folder_x<min_tt, percol_, simd_<avx_t>, bloc, M, N>();
}

MN_CASE( min_percol_avx_grid )
{
	test_folder_x<min_tt, percol_, simd_<avx_t>, grid, M, N>();
}

#endif

MN_CASE( min_auto_cont )
{
	test_folder<min_tt, cont, M, N>();
}

MN_CASE( min_auto_bloc )
{
	test_folder<min_tt, bloc, M, N>();
}

MN_CASE( min_auto_grid )
{
	test_folder<min_tt, grid, M, N>();
}



// Packs

// sum

AUTO_TPACK( sum_linear_scalar )
{
	ADD_MN_CASE_3X3( sum_linear_scalar_cont, DM, DN )

	ADD_MN_CASE( sum_linear_scalar_bloc, 1, 1 )
	ADD_MN_CASE( sum_linear_scalar_bloc, 1, 0 )
	ADD_MN_CASE( sum_linear_scalar_bloc, 1, DN )
	ADD_MN_CASE( sum_linear_scalar_bloc, 0, 1 )
	ADD_MN_CASE( sum_linear_scalar_bloc, DM, 1 )

	ADD_MN_CASE( sum_linear_scalar_grid, 1, 1 )
	ADD_MN_CASE( sum_linear_scalar_grid, 1, 0 )
	ADD_MN_CASE( sum_linear_scalar_grid, 1, DN )
	ADD_MN_CASE( sum_linear_scalar_grid, 0, 1 )
	ADD_MN_CASE( sum_linear_scalar_grid, DM, 1 )
}

AUTO_TPACK( sum_linear_sse )
{
	ADD_MN_CASE_3X3( sum_linear_sse_cont, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( sum_linear_avx )
{
	ADD_MN_CASE_3X3( sum_linear_avx_cont, DM, DN )
}
#endif

AUTO_TPACK( sum_percol_scalar )
{
	ADD_MN_CASE_3X3( sum_percol_scalar_cont, DM, DN )
	ADD_MN_CASE_3X3( sum_percol_scalar_bloc, DM, DN )
	ADD_MN_CASE_3X3( sum_percol_scalar_grid, DM, DN )
}

AUTO_TPACK( sum_percol_sse )
{
	ADD_MN_CASE_3X3( sum_percol_sse_cont, DM, DN )
	ADD_MN_CASE_3X3( sum_percol_sse_bloc, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( sum_percol_avx )
{
	ADD_MN_CASE_3X3( sum_percol_avx_cont, DM, DN )
	ADD_MN_CASE_3X3( sum_percol_avx_bloc, DM, DN )
}
#endif

AUTO_TPACK( sum_auto )
{
	ADD_MN_CASE_3X3( sum_auto_cont, DM, DN )
	ADD_MN_CASE_3X3( sum_auto_bloc, DM, DN )
	ADD_MN_CASE_3X3( sum_auto_grid, DM, DN )
}


// max

AUTO_TPACK( max_linear_scalar )
{
	ADD_MN_CASE_3X3( max_linear_scalar_cont, DM, DN )

	ADD_MN_CASE( max_linear_scalar_bloc, 1, 1 )
	ADD_MN_CASE( max_linear_scalar_bloc, 1, 0 )
	ADD_MN_CASE( max_linear_scalar_bloc, 1, DN )
	ADD_MN_CASE( max_linear_scalar_bloc, 0, 1 )
	ADD_MN_CASE( max_linear_scalar_bloc, DM, 1 )

	ADD_MN_CASE( max_linear_scalar_grid, 1, 1 )
	ADD_MN_CASE( max_linear_scalar_grid, 1, 0 )
	ADD_MN_CASE( max_linear_scalar_grid, 1, DN )
	ADD_MN_CASE( max_linear_scalar_grid, 0, 1 )
	ADD_MN_CASE( max_linear_scalar_grid, DM, 1 )
}

AUTO_TPACK( max_linear_sse )
{
	ADD_MN_CASE_3X3( max_linear_sse_cont, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( max_linear_avx )
{
	ADD_MN_CASE_3X3( max_linear_avx_cont, DM, DN )
}
#endif

AUTO_TPACK( max_percol_scalar )
{
	ADD_MN_CASE_3X3( max_percol_scalar_cont, DM, DN )
	ADD_MN_CASE_3X3( max_percol_scalar_bloc, DM, DN )
	ADD_MN_CASE_3X3( max_percol_scalar_grid, DM, DN )
}

AUTO_TPACK( max_percol_sse )
{
	ADD_MN_CASE_3X3( max_percol_sse_cont, DM, DN )
	ADD_MN_CASE_3X3( max_percol_sse_bloc, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( max_percol_avx )
{
	ADD_MN_CASE_3X3( max_percol_avx_cont, DM, DN )
	ADD_MN_CASE_3X3( max_percol_avx_bloc, DM, DN )
}
#endif

AUTO_TPACK( max_auto )
{
	ADD_MN_CASE_3X3( max_auto_cont, DM, DN )
	ADD_MN_CASE_3X3( max_auto_bloc, DM, DN )
	ADD_MN_CASE_3X3( max_auto_grid, DM, DN )
}


// min

AUTO_TPACK( min_linear_scalar )
{
	ADD_MN_CASE_3X3( min_linear_scalar_cont, DM, DN )

	ADD_MN_CASE( min_linear_scalar_bloc, 1, 1 )
	ADD_MN_CASE( min_linear_scalar_bloc, 1, 0 )
	ADD_MN_CASE( min_linear_scalar_bloc, 1, DN )
	ADD_MN_CASE( min_linear_scalar_bloc, 0, 1 )
	ADD_MN_CASE( min_linear_scalar_bloc, DM, 1 )

	ADD_MN_CASE( min_linear_scalar_grid, 1, 1 )
	ADD_MN_CASE( min_linear_scalar_grid, 1, 0 )
	ADD_MN_CASE( min_linear_scalar_grid, 1, DN )
	ADD_MN_CASE( min_linear_scalar_grid, 0, 1 )
	ADD_MN_CASE( min_linear_scalar_grid, DM, 1 )
}

AUTO_TPACK( min_linear_sse )
{
	ADD_MN_CASE_3X3( min_linear_sse_cont, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( min_linear_avx )
{
	ADD_MN_CASE_3X3( min_linear_avx_cont, DM, DN )
}
#endif

AUTO_TPACK( min_percol_scalar )
{
	ADD_MN_CASE_3X3( min_percol_scalar_cont, DM, DN )
	ADD_MN_CASE_3X3( min_percol_scalar_bloc, DM, DN )
	ADD_MN_CASE_3X3( min_percol_scalar_grid, DM, DN )
}

AUTO_TPACK( min_percol_sse )
{
	ADD_MN_CASE_3X3( min_percol_sse_cont, DM, DN )
	ADD_MN_CASE_3X3( min_percol_sse_bloc, DM, DN )
}

#ifdef LMAT_HAS_AVX
AUTO_TPACK( min_percol_avx )
{
	ADD_MN_CASE_3X3( min_percol_avx_cont, DM, DN )
	ADD_MN_CASE_3X3( min_percol_avx_bloc, DM, DN )
}
#endif

AUTO_TPACK( min_auto )
{
	ADD_MN_CASE_3X3( min_auto_cont, DM, DN )
	ADD_MN_CASE_3X3( min_auto_bloc, DM, DN )
	ADD_MN_CASE_3X3( min_auto_grid, DM, DN )
}


