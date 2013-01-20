/**
 * @file test_sample_wor.cpp
 *
 * @brief Unit testing of sampling without replacement
 *
 * @author Dahua Lin
 */


#include "distr_test_base.h"
#include <light_mat/random/sample_wor.h>


default_rand_stream rstream;
const index_t N = 10000;
const index_t L = 5;

template<class Enumerator>
void test_sample_wor(Enumerator& enumerator)
{
	typedef typename Enumerator::value_type T;

	ASSERT_EQ( enumerator.length(), L );

	dense_matrix<uint32_t> cnts(L, L, zero());
	dense_col<bool> visited(L);

	for (index_t i = 0; i < N; ++i)
	{
		enumerator.reset();
		ASSERT_EQ( enumerator.remain(), L );

		zero(visited);

		for (index_t j = 0; j < L; ++j)
		{
			ASSERT_FALSE( enumerator.is_end() );

			T x = enumerator.next();

			index_t ix = (index_t)x;
			ASSERT_TRUE( ix >= 0 && ix < L );

			bool x_visited = visited[ix];
			ASSERT_FALSE( x_visited );
			visited[ix] = true;

			++ cnts(j, ix);

			ASSERT_EQ( enumerator.remain(), L-j-1 );
		}

		ASSERT_TRUE( enumerator.is_end() );
	}

	dense_matrix<double> p0(L, L, zero());
	fill(p0, 1.0 / double(L));

	dense_matrix<double> p(L, L);
	for (index_t j = 0; j < L; ++j)
	{
		for (index_t i = 0; i < L; ++i)
			p(i, j) = double(cnts(i, j)) / double(N);
	}

	double tol = get_p_tol(N);
	ASSERT_MAT_APPROX(L, L, p, p0, tol);
}


T_CASE( test_shuffler )
{
	rand_shuffle_enumerator<T> enumerator(rstream, L);
	test_sample_wor(enumerator);
}

T_CASE( test_past_avoider )
{
	past_avoid_rand_enumerator<T> enumerator(rstream, L);
	test_sample_wor(enumerator);
}


AUTO_TPACK( test_shuffler )
{
	ADD_T_CASE( test_shuffler, uint32_t )
	ADD_T_CASE( test_shuffler, int32_t )
}


AUTO_TPACK( test_past_avoider )
{
	ADD_T_CASE( test_past_avoider, uint32_t )
	ADD_T_CASE( test_past_avoider, int32_t )
}

