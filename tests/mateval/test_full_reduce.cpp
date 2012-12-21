/**
 * @file test_full_reduce.cpp
 *
 * Test of full reduction on matrices
 * 
 * @author Dahua Lin 
 */

#include "../test_base.h"
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/mat_reduce.h>
#include <cstdlib>

using namespace lmat;
using namespace lmat::test;

inline double randunif()
{
	double u = (double)std::rand() / double(RAND_MAX);
	return u * 2.0 - 1.0;
}

template<class Mat, typename T>
void fill_rand(IRegularMatrix<Mat, T>& mat)
{
	for (index_t j = 0; j < mat.ncolumns(); ++j)
	{
		for (index_t i = 0; i < mat.nrows(); ++i)
		{
			mat(i, j) = randunif();
		}
	}
}

const index_t max_len = 48;

SIMPLE_CASE( full_reduce, sum )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += s[i];

		double r = sum(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, mean )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += s[i];
		r0 /= k;

		double r = mean(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


SIMPLE_CASE( full_reduce, maximum )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		double r0 = - std::numeric_limits<double>::infinity();
		for (index_t i = 0; i < k; ++i) r0 = math::max(r0, s[i]);

		double r = maximum(s(range(0, k)));
		ASSERT_EQ(r, r0);
	}
}

SIMPLE_CASE( full_reduce, minimum )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 1; k <= max_len; ++k)
	{
		double r0 = std::numeric_limits<double>::infinity();
		for (index_t i = 0; i < k; ++i) r0 = math::min(r0, s[i]);

		double r = minimum(s(range(0, k)));
		ASSERT_EQ(r, r0);
	}
}


SIMPLE_CASE( full_reduce, asum )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(s[i]);

		double r = asum(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, amean )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(s[i]);
		r0 /= k;

		double r = amean(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


SIMPLE_CASE( full_reduce, amax )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i)
			r0 = math::max(r0, math::abs(s[i]));

		double r = amax(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, sqsum )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::sqr(s[i]);

		double r = sqsum(s(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


SIMPLE_CASE( full_reduce, diff_asum )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(s[i] - s2[i]);

		double r = diff_asum(s(range(0, k)), s2(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, diff_amean )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 1; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::abs(s[i] - s2[i]);
		r0 /= k;

		double r = diff_amean(s(range(0, k)), s2(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, diff_amax )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i)
			r0 = math::max(r0, math::abs(s[i] - s2[i]));

		double r = diff_amax(s(range(0, k)), s2(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}

SIMPLE_CASE( full_reduce, diff_sqsum )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += math::sqr(s[i] - s2[i]);

		double r = diff_sqsum(s(range(0, k)), s2(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


SIMPLE_CASE( full_reduce, dot )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 0; k <= max_len; ++k)
	{
		double r0 = 0;
		for (index_t i = 0; i < k; ++i) r0 += s[i] * s2[i];

		double r = dot(s(range(0, k)), s2(range(0, k)));
		ASSERT_APPROX(r, r0, 1.0e-12);
	}
}


SIMPLE_CASE( full_reduce, norms )
{
	dense_col<double> s(max_len);
	fill_rand(s);

	for (index_t k = 0; k <= max_len; ++k)
	{
		ref_col<double> sk = s(range(0, k));

		ASSERT_EQ( norm(sk, norms::L1_()), asum(sk) );
		ASSERT_EQ( norm(sk, norms::L2_()), math::sqrt(sqsum(sk)) );
		ASSERT_EQ( norm(sk, norms::Linf_()), amax(sk) );
	}
}


SIMPLE_CASE( full_reduce, diff_norms )
{
	dense_col<double> s(max_len);
	dense_col<double> s2(max_len);
	fill_rand(s);
	fill_rand(s2);

	for (index_t k = 0; k <= max_len; ++k)
	{
		ref_col<double> sk = s(range(0, k));
		ref_col<double> sk2 = s2(range(0, k));

		ASSERT_EQ( diff_norm(sk, sk2, norms::L1_()), diff_asum(sk, sk2) );
		ASSERT_EQ( diff_norm(sk, sk2, norms::L2_()), math::sqrt(diff_sqsum(sk, sk2)) );
		ASSERT_EQ( diff_norm(sk, sk2, norms::Linf_()), diff_amax(sk, sk2) );
	}
}





BEGIN_TPACK( full_reduce )
	ADD_SIMPLE_CASE( full_reduce, sum )
	ADD_SIMPLE_CASE( full_reduce, mean )
	ADD_SIMPLE_CASE( full_reduce, maximum )
	ADD_SIMPLE_CASE( full_reduce, minimum )

	ADD_SIMPLE_CASE( full_reduce, asum )
	ADD_SIMPLE_CASE( full_reduce, amean )
	ADD_SIMPLE_CASE( full_reduce, amax )
	ADD_SIMPLE_CASE( full_reduce, sqsum )

	ADD_SIMPLE_CASE( full_reduce, diff_asum )
	ADD_SIMPLE_CASE( full_reduce, diff_amean )
	ADD_SIMPLE_CASE( full_reduce, diff_amax )
	ADD_SIMPLE_CASE( full_reduce, diff_sqsum )

	ADD_SIMPLE_CASE( full_reduce, dot )

	ADD_SIMPLE_CASE( full_reduce, norms )
	ADD_SIMPLE_CASE( full_reduce, diff_norms )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( full_reduce )
END_MAIN_SUITE

