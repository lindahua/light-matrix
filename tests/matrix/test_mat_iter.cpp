/**
 * @file test_mat_iter.cpp
 *
 * @brief Unit testing of matrix iteration
 *
 * @author Dahua Lin
 */

#include "../test_base.h"
#include "../multimat_supp.h"

#include <light_mat/matrix/matrix_classes.h>
#include <vector>

using namespace lmat;
using namespace lmat::test;


template<class MTag, int M, int N>
void test_matrix_matiter()
{
	typedef mat_host<MTag, double, M, N> host_t;

	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	host_t host(m, n);
	host.fill_lin();

	mat_t a = host.get_mat();
	cmat_t c = host.get_cmat();
	dense_matrix<double> r(a);

	std::vector<double> av(begin(a), end(a));
	ASSERT_EQ( (index_t)av.size(), m * n );
	ASSERT_VEC_EQ( m * n, &av[0], r );

	std::vector<double> cv(begin(c), end(c));
	ASSERT_EQ( (index_t)cv.size(), m * n );
	ASSERT_VEC_EQ( m * n, &cv[0], r );

	double v0 = 2.5;
	*begin(a) = v0;
	ASSERT_EQ(a(0, 0), v0);
}


template<class MTag, int M, int N>
void test_matrix_coliter()
{
	typedef mat_host<MTag, double, M, N> host_t;

	typedef typename host_t::cmat_t cmat_t;
	typedef typename host_t::mat_t mat_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	host_t host(m, n);
	host.fill_lin();

	mat_t a = host.get_mat();
	cmat_t c = host.get_cmat();
	dense_matrix<double> r(a);

	for (index_t j = 0; j < n; ++j)
	{
		auto rj = a.column(j);

		std::vector<double> av(a.col_begin(j), a.col_end(j));
		ASSERT_EQ( (index_t)av.size(), m );
		ASSERT_VEC_EQ( m, &av[0], rj );

		std::vector<double> cv(c.col_begin(j), c.col_end(j));
		ASSERT_EQ( (index_t)cv.size(), m );
		ASSERT_VEC_EQ( m, &cv[0], rj );
	}
}


MN_CASE( mat_iter, cont )
{
	test_matrix_matiter<cont, M, N>();
}

MN_CASE( mat_iter, bloc )
{
	test_matrix_matiter<bloc, M, N>();
}

MN_CASE( mat_iter, grid )
{
	test_matrix_matiter<grid, M, N>();
}

MN_CASE( mat_coliter, cont )
{
	test_matrix_coliter<cont, M, N>();
}

MN_CASE( mat_coliter, bloc )
{
	test_matrix_coliter<bloc, M, N>();
}

MN_CASE( mat_coliter, grid )
{
	test_matrix_coliter<grid, M, N>();
}


BEGIN_TPACK( mat_iter_cont )
	ADD_MN_CASE_3X3( mat_iter, cont, DM, DN )
END_TPACK

BEGIN_TPACK( mat_iter_bloc )
	ADD_MN_CASE_3X3( mat_iter, bloc, DM, DN )
END_TPACK

BEGIN_TPACK( mat_iter_grid )
	ADD_MN_CASE_3X3( mat_iter, grid, DM, DN )
END_TPACK

BEGIN_TPACK( mat_coliter_cont )
	ADD_MN_CASE_3X3( mat_coliter, cont, DM, DN )
END_TPACK

BEGIN_TPACK( mat_coliter_bloc )
	ADD_MN_CASE_3X3( mat_coliter, bloc, DM, DN )
END_TPACK

BEGIN_TPACK( mat_coliter_grid )
	ADD_MN_CASE_3X3( mat_coliter, grid, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( mat_iter_cont )
	ADD_TPACK( mat_iter_bloc )
	ADD_TPACK( mat_iter_grid )
	ADD_TPACK( mat_coliter_cont )
	ADD_TPACK( mat_coliter_bloc )
	ADD_TPACK( mat_coliter_grid )
END_MAIN_SUITE

