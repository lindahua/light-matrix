/**
 * @file test_matrix_props.cpp
 *
 * @brief Test matrix property extraction
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include <light_mat/matrix/dense_matrix.h>
#include <light_mat/matrix/ref_matrix.h>
#include <light_mat/matrix/ref_block.h>
#include <light_mat/matrix/ref_grid.h>

using namespace lmat;

MN_CASE( mat_ct_dims )
{
	typedef dense_matrix<double, M, N> mat_t;

	ASSERT_CT_VALUE( meta::nrows<mat_t>, M );
	ASSERT_CT_VALUE( meta::ncols<mat_t>, N );
	ASSERT_CT_VALUE( meta::nelems<mat_t>, M * N );

	typedef typename meta::shape<mat_t>::type shape_t;
	typedef matrix_shape<M, N> shape_t_expect;
	ASSERT_SAME_TYPE( shape_t, shape_t_expect );
}


MN_CASE( mat_ct_layout )
{
	typedef dense_matrix<double, M, N> dens_t;
	typedef ref_matrix<double, M, N> cont_t;
	typedef ref_block<double, M, N> bloc_t;
	typedef ref_grid<double, M, N> grid_t;

	bool cont_is_continuous = true;
	bool bloc_is_continuous = N == 1;
	bool grid_is_continuous = M == 1 && N == 1;

	bool cont_is_pc_continuous = true;
	bool bloc_is_pc_continuous = true;
	bool grid_is_pc_continuous = M == 1;

	ASSERT_CT_VALUE( meta::is_contiguous<dens_t>, cont_is_continuous );
	ASSERT_CT_VALUE( meta::is_contiguous<cont_t>, cont_is_continuous );
	ASSERT_CT_VALUE( meta::is_contiguous<bloc_t>, bloc_is_continuous );
	ASSERT_CT_VALUE( meta::is_contiguous<grid_t>, grid_is_continuous );

	ASSERT_CT_VALUE( meta::is_percol_contiguous<dens_t>, cont_is_pc_continuous );
	ASSERT_CT_VALUE( meta::is_percol_contiguous<cont_t>, cont_is_pc_continuous );
	ASSERT_CT_VALUE( meta::is_percol_contiguous<bloc_t>, bloc_is_pc_continuous );
	ASSERT_CT_VALUE( meta::is_percol_contiguous<grid_t>, grid_is_pc_continuous );
}


// compatible rows

template<int M, int N>
struct binary_compatible_nrows
{
	typedef dense_matrix<double, M, 1> A1;
	typedef dense_matrix<double, N, 1> A2;

	static const bool value = meta::have_compatible_nrows<A1, A2>::value;
};

template<int M, int N, int K>
struct ternary_compatible_nrows
{
	typedef dense_matrix<double, M, 1> A1;
	typedef dense_matrix<double, N, 1> A2;
	typedef dense_matrix<double, K, 1> A3;

	static const bool value = meta::have_compatible_nrows<A1, A2, A3>::value;
};


template<int M>
struct unary_common_nrows
{
	typedef dense_matrix<double, M, 1> A1;

	static const int value = meta::common_nrows<A1>::value;
};

template<int M, int N>
struct binary_common_nrows
{
	typedef dense_matrix<double, M, 1> A1;
	typedef dense_matrix<double, N, 1> A2;

	static const int value = meta::common_nrows<A1, A2>::value;
};

template<int M, int N, int K>
struct ternary_common_nrows
{
	typedef dense_matrix<double, M, 1> A1;
	typedef dense_matrix<double, N, 1> A2;
	typedef dense_matrix<double, K, 1> A3;

	static const int value = meta::common_nrows<A1, A2, A3>::value;
};


SIMPLE_CASE( mat_ct_compat_nrows )
{
	typedef binary_compatible_nrows<0, 3> t00;
	ASSERT_CT_VALUE( t00,  true );

	typedef binary_compatible_nrows<0, 3> t03;
	ASSERT_CT_VALUE( t03,  true );

	typedef binary_compatible_nrows<3, 0> t30;
	ASSERT_CT_VALUE( t30,  true );

	typedef binary_compatible_nrows<3, 3> t33;
	ASSERT_CT_VALUE( t33,  true );

	typedef binary_compatible_nrows<3, 4> t34;
	ASSERT_CT_VALUE( t34,  false );

	typedef ternary_compatible_nrows<0, 0, 0> t000;
	ASSERT_CT_VALUE( t000,  true );

	typedef ternary_compatible_nrows<0, 0, 3> t003;
	ASSERT_CT_VALUE( t003,  true );

	typedef ternary_compatible_nrows<0, 3, 0> t030;
	ASSERT_CT_VALUE( t030,  true );

	typedef ternary_compatible_nrows<0, 3, 3> t033;
	ASSERT_CT_VALUE( t033,  true );

	typedef ternary_compatible_nrows<3, 0, 0> t300;
	ASSERT_CT_VALUE( t300,  true );

	typedef ternary_compatible_nrows<3, 0, 3> t303;
	ASSERT_CT_VALUE( t303,  true );

	typedef ternary_compatible_nrows<3, 3, 0> t330;
	ASSERT_CT_VALUE( t330,  true );

	typedef ternary_compatible_nrows<3, 3, 3> t333;
	ASSERT_CT_VALUE( t333,  true );

	typedef ternary_compatible_nrows<0, 3, 4> t034;
	ASSERT_CT_VALUE( t034,  false );

	typedef ternary_compatible_nrows<3, 4, 0> t340;
	ASSERT_CT_VALUE( t340,  false );

	typedef ternary_compatible_nrows<3, 4, 5> t345;
	ASSERT_CT_VALUE( t345,  false );
}


SIMPLE_CASE( mat_ct_common_nrows )
{
	typedef unary_common_nrows<0> t0;
	ASSERT_CT_VALUE( t0, 0 );

	typedef unary_common_nrows<3> t3;
	ASSERT_CT_VALUE( t3, 3 );

	typedef binary_common_nrows<0, 0> t00;
	ASSERT_CT_VALUE( t00,  0 );

	typedef binary_common_nrows<0, 3> t03;
	ASSERT_CT_VALUE( t03,  3 );

	typedef binary_common_nrows<3, 0> t30;
	ASSERT_CT_VALUE( t30,  3 );

	typedef binary_common_nrows<3, 3> t33;
	ASSERT_CT_VALUE( t33,  3 );

	typedef ternary_common_nrows<0, 0, 0> t000;
	ASSERT_CT_VALUE( t000,  0 );

	typedef ternary_common_nrows<0, 0, 3> t003;
	ASSERT_CT_VALUE( t003,  3 );

	typedef ternary_common_nrows<0, 3, 0> t030;
	ASSERT_CT_VALUE( t030,  3 );

	typedef ternary_common_nrows<0, 3, 3> t033;
	ASSERT_CT_VALUE( t033,  3 );

	typedef ternary_common_nrows<3, 0, 0> t300;
	ASSERT_CT_VALUE( t300,  3 );

	typedef ternary_common_nrows<3, 0, 3> t303;
	ASSERT_CT_VALUE( t303,  3 );

	typedef ternary_common_nrows<3, 3, 0> t330;
	ASSERT_CT_VALUE( t330,  3 );

	typedef ternary_common_nrows<3, 3, 3> t333;
	ASSERT_CT_VALUE( t333,  3 );
}


// compatible cols

template<int M, int N>
struct binary_compatible_ncols
{
	typedef dense_matrix<double, 1, M> A1;
	typedef dense_matrix<double, 1, N> A2;

	static const bool value = meta::have_compatible_ncols<A1, A2>::value;
};

template<int M, int N, int K>
struct ternary_compatible_ncols
{
	typedef dense_matrix<double, 1, M> A1;
	typedef dense_matrix<double, 1, N> A2;
	typedef dense_matrix<double, 1, K> A3;

	static const bool value = meta::have_compatible_ncols<A1, A2, A3>::value;
};


template<int M>
struct unary_common_ncols
{
	typedef dense_matrix<double, 1, M> A1;

	static const int value = meta::common_ncols<A1>::value;
};

template<int M, int N>
struct binary_common_ncols
{
	typedef dense_matrix<double, 1, M> A1;
	typedef dense_matrix<double, 1, N> A2;

	static const int value = meta::common_ncols<A1, A2>::value;
};

template<int M, int N, int K>
struct ternary_common_ncols
{
	typedef dense_matrix<double, 1, M> A1;
	typedef dense_matrix<double, 1, N> A2;
	typedef dense_matrix<double, 1, K> A3;

	static const int value = meta::common_ncols<A1, A2, A3>::value;
};



SIMPLE_CASE( mat_ct_compat_ncols )
{
	typedef binary_compatible_ncols<0, 3> t00;
	ASSERT_CT_VALUE( t00,  true );

	typedef binary_compatible_ncols<0, 3> t03;
	ASSERT_CT_VALUE( t03,  true );

	typedef binary_compatible_ncols<3, 0> t30;
	ASSERT_CT_VALUE( t30,  true );

	typedef binary_compatible_ncols<3, 3> t33;
	ASSERT_CT_VALUE( t33,  true );

	typedef binary_compatible_ncols<3, 4> t34;
	ASSERT_CT_VALUE( t34,  false );

	typedef ternary_compatible_ncols<0, 0, 0> t000;
	ASSERT_CT_VALUE( t000,  true );

	typedef ternary_compatible_ncols<0, 0, 3> t003;
	ASSERT_CT_VALUE( t003,  true );

	typedef ternary_compatible_ncols<0, 3, 0> t030;
	ASSERT_CT_VALUE( t030,  true );

	typedef ternary_compatible_ncols<0, 3, 3> t033;
	ASSERT_CT_VALUE( t033,  true );

	typedef ternary_compatible_ncols<3, 0, 0> t300;
	ASSERT_CT_VALUE( t300,  true );

	typedef ternary_compatible_ncols<3, 0, 3> t303;
	ASSERT_CT_VALUE( t303,  true );

	typedef ternary_compatible_ncols<3, 3, 0> t330;
	ASSERT_CT_VALUE( t330,  true );

	typedef ternary_compatible_ncols<3, 3, 3> t333;
	ASSERT_CT_VALUE( t333,  true );

	typedef ternary_compatible_ncols<0, 3, 4> t034;
	ASSERT_CT_VALUE( t034,  false );

	typedef ternary_compatible_ncols<3, 4, 0> t340;
	ASSERT_CT_VALUE( t340,  false );

	typedef ternary_compatible_ncols<3, 4, 5> t345;
	ASSERT_CT_VALUE( t345,  false );
}


SIMPLE_CASE( mat_ct_common_ncols )
{
	typedef unary_common_ncols<0> t0;
	ASSERT_CT_VALUE( t0, 0 );

	typedef unary_common_ncols<3> t3;
	ASSERT_CT_VALUE( t3, 3 );

	typedef binary_common_ncols<0, 0> t00;
	ASSERT_CT_VALUE( t00,  0 );

	typedef binary_common_ncols<0, 3> t03;
	ASSERT_CT_VALUE( t03,  3 );

	typedef binary_common_ncols<3, 0> t30;
	ASSERT_CT_VALUE( t30,  3 );

	typedef binary_common_ncols<3, 3> t33;
	ASSERT_CT_VALUE( t33,  3 );

	typedef ternary_common_ncols<0, 0, 0> t000;
	ASSERT_CT_VALUE( t000,  0 );

	typedef ternary_common_ncols<0, 0, 3> t003;
	ASSERT_CT_VALUE( t003,  3 );

	typedef ternary_common_ncols<0, 3, 0> t030;
	ASSERT_CT_VALUE( t030,  3 );

	typedef ternary_common_ncols<0, 3, 3> t033;
	ASSERT_CT_VALUE( t033,  3 );

	typedef ternary_common_ncols<3, 0, 0> t300;
	ASSERT_CT_VALUE( t300,  3 );

	typedef ternary_common_ncols<3, 0, 3> t303;
	ASSERT_CT_VALUE( t303,  3 );

	typedef ternary_common_ncols<3, 3, 0> t330;
	ASSERT_CT_VALUE( t330,  3 );

	typedef ternary_common_ncols<3, 3, 3> t333;
	ASSERT_CT_VALUE( t333,  3 );
}

LTEST_INIT_AUTOSUITE

AUTO_TPACK( mat_ct_dims )
{
	ADD_MN_CASE_3X3( mat_ct_dims, 3, 4 )
}

AUTO_TPACK( mat_ct_layout )
{
	ADD_MN_CASE_3X3( mat_ct_layout, 3, 4 )
}

AUTO_TPACK( mat_common_dims )
{
	ADD_SIMPLE_CASE( mat_ct_compat_nrows )
	ADD_SIMPLE_CASE( mat_ct_common_nrows )
	ADD_SIMPLE_CASE( mat_ct_compat_ncols )
	ADD_SIMPLE_CASE( mat_ct_common_ncols )
}




