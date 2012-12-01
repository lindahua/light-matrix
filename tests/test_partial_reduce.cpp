/**
 * @file test_partial_reduce.cpp
 *
 * Unit testing of partial reduction
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/partial_reduce.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 6;
const index_t DN = 7;
const index_t LDim = 8;

MN_CASE( colwise_reduce, sum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j);
		r0[j] = s;
	}

	// test

	row_t r = sum(A, colwise());

	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );

}

MN_CASE( rowwise_reduce, sum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j);
		r0[i] = s;
	}

	// test

	col_t r = sum(A, rowwise());

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );

	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, mean )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j);
		r0[j] = s / double(m);
	}

	// test

	row_t r = mean(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_APPROX( n, r, r0, 1.0e-12 );
}


MN_CASE( rowwise_reduce, mean )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j);
		r0[i] = s / double(n);
	}

	// test

	col_t r = mean(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_APPROX( m, r, r0, 1.0e-12 );
}


MN_CASE( colwise_reduce, maximum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (10 - i));

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = A(0, j);
		for (index_t i = 0; i < m; ++i)
			if (A(i, j) > s) s = A(i, j);
		r0[j] = s;
	}

	// test

	row_t r = maximum(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, maximum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (10 - i));

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = A(i, 0);
		for (index_t j = 0; j < n; ++j)
			if (A(i, j) > s) s = A(i, j);
		r0[i] = s;
	}

	// test

	col_t r = maximum(A, rowwise());

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, minimum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (i - 10));

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = A(0, j);
		for (index_t i = 0; i < m; ++i)
			if (A(i, j) < s) s = A(i, j);
		r0[j] = s;
	}

	// test

	row_t r = minimum(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, minimum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double((i + 1) * (i - 10));

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = A(i, 0);
		for (index_t j = 0; j < n; ++j)
			if (A(i, j) < s) s = A(i, j);
		r0[i] = s;
	}

	// test

	col_t r = minimum(A, rowwise());

	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}





MN_CASE( colwise_reduce, L1norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = i & 1 ? double(i + 1) : double(-(i + 1));
	}

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += std::abs(A(i, j));
		r0[j] = s;
	}

	// test

	row_t r = L1norm(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, L1norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = i & 1 ? double(i + 1) : double(-(i + 1));
	}

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += std::abs(A(i, j));
		r0[i] = s;
	}

	// test

	col_t r = L1norm(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}



MN_CASE( colwise_reduce, sqL2norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j) * A(i, j);
		r0[j] = s;
	}

	// test

	row_t r = sqL2norm(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, sqL2norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j) * A(i, j);
		r0[i] = s;
	}

	// test

	col_t r = sqL2norm(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, L2norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j) * A(i, j);
		r0[j] = std::sqrt(s);
	}

	// test

	row_t r = L2norm(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_APPROX( n, r, r0, 1.0e-12 );
}

MN_CASE( rowwise_reduce, L2norm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j) * A(i, j);
		r0[i] = std::sqrt(s);
	}

	// test

	col_t r = L2norm(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_APPROX( m, r, r0, 1.0e-12 );
}


MN_CASE( colwise_reduce, Linfnorm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = i & 1 ? double(i + 1) : double(-(i + 1));
	}

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i)
		{
			if ( std::abs(A(i, j)) > s ) s = std::abs(A(i, j));
		}
		r0[j] = s;
	}

	// test

	row_t r = Linfnorm(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, Linfnorm )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = i & 1 ? double(i + 1) : double(-(i + 1));
	}

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j)
		{
			if ( std::abs(A(i, j)) > s ) s = std::abs(A(i, j));
		}
		r0[i] = s;
	}

	// test

	col_t r = Linfnorm(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, logsum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 2);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += std::log(A(i, j));
		r0[j] = s;
	}

	// test

	row_t r = logsum(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, logsum )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += std::log(A(i, j));
		r0[i] = s;
	}

	// test

	col_t r = logsum(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, entropy )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 2);

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s -= A(i, j) * std::log(A(i, j));
		r0[j] = s;
	}

	// test

	row_t r = entropy(A, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}

MN_CASE( rowwise_reduce, entropy )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	for (index_t i = 0; i < m * n; ++i) A[i] = double(i + 1);

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s -= A(i, j) * std::log(A(i, j));
		r0[i] = s;
	}

	// test

	col_t r = entropy(A, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}



MN_CASE( colwise_reduce, dot )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = double(i + 1);
		B[i] = double(i + 2);
	}

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double s = 0.0;
		for (index_t i = 0; i < m; ++i) s += A(i, j) * B(i, j);
		r0[j] = s;
	}

	// test

	row_t r = dot(A, B, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_EQ( n, r, r0 );
}


MN_CASE( rowwise_reduce, dot )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = double(i + 1);
		B[i] = double(i + 2);
	}

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double s = 0.0;
		for (index_t j = 0; j < n; ++j) s += A(i, j) * B(i, j);
		r0[i] = s;
	}

	// test

	col_t r = dot(A, B, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_EQ( m, r, r0 );
}


MN_CASE( colwise_reduce, nrmdot )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_row<double, N> row_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = double(i + 1);
		B[i] = double(2 * i + 1);
	}

	// prepare ground-truth

	row_t r0(n);
	for (index_t j = 0; j < n; ++j)
	{
		double sab = 0.0;
		double saa = 0.0;
		double sbb = 0.0;

		for (index_t i = 0; i < m; ++i)
		{
			sab += A(i, j) * B(i, j);
			saa += A(i, j) * A(i, j);
			sbb += B(i, j) * B(i, j);
		}
		r0[j] = sab / (sqrt(saa) * sqrt(sbb));
	}

	// test

	row_t r = nrmdot(A, B, colwise());
	ASSERT_EQ( r.nrows(), 1 );
	ASSERT_EQ( r.ncolumns(), n );
	ASSERT_VEC_APPROX( n, r, r0, 1.0e-12 );
}

MN_CASE( rowwise_reduce, nrmdot )
{
	typedef dense_matrix<double, M, N> mat_t;
	typedef dense_col<double, M> col_t;

	const index_t m = M == 0 ? DM : M;
	const index_t n = N == 0 ? DN : N;

	mat_t A(m, n);
	mat_t B(m, n);
	for (index_t i = 0; i < m * n; ++i)
	{
		A[i] = double(i + 1);
		B[i] = double(2 * i + 1);
	}

	// prepare ground-truth

	col_t r0(m);
	for (index_t i = 0; i < m; ++i)
	{
		double sab = 0.0;
		double saa = 0.0;
		double sbb = 0.0;

		for (index_t j = 0; j < n; ++j)
		{
			sab += A(i, j) * B(i, j);
			saa += A(i, j) * A(i, j);
			sbb += B(i, j) * B(i, j);
		}
		r0[i] = sab / (sqrt(saa) * sqrt(sbb));
	}

	// test

	col_t r = nrmdot(A, B, rowwise());
	ASSERT_EQ( r.nrows(), m );
	ASSERT_EQ( r.ncolumns(), 1 );
	ASSERT_VEC_APPROX( m, r, r0, 1.0e-12 );
}


BEGIN_TPACK( colwise_sum )
	ADD_MN_CASE_3X3( colwise_reduce, sum, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_sum )
	ADD_MN_CASE_3X3( rowwise_reduce, sum, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_mean )
	ADD_MN_CASE_3X3( colwise_reduce, mean, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_mean )
	ADD_MN_CASE_3X3( rowwise_reduce, mean, DM, DN )
END_TPACK


BEGIN_TPACK( colwise_maximum )
	ADD_MN_CASE_3X3( colwise_reduce, maximum, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_maximum )
	ADD_MN_CASE_3X3( rowwise_reduce, maximum, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_minimum )
	ADD_MN_CASE_3X3( colwise_reduce, minimum, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_minimum )
	ADD_MN_CASE_3X3( rowwise_reduce, minimum, DM, DN )
END_TPACK


BEGIN_TPACK( colwise_L1norm )
	ADD_MN_CASE_3X3( colwise_reduce, L1norm, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_L1norm )
	ADD_MN_CASE_3X3( rowwise_reduce, L1norm, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_sqL2norm )
	ADD_MN_CASE_3X3( colwise_reduce, sqL2norm, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_sqL2norm )
	ADD_MN_CASE_3X3( rowwise_reduce, sqL2norm, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_L2norm )
	ADD_MN_CASE_3X3( colwise_reduce, L2norm, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_L2norm )
	ADD_MN_CASE_3X3( rowwise_reduce, L2norm, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_Linfnorm )
	ADD_MN_CASE_3X3( colwise_reduce, Linfnorm, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_Linfnorm )
	ADD_MN_CASE_3X3( rowwise_reduce, Linfnorm, DM, DN )
END_TPACK


BEGIN_TPACK( colwise_logsum )
	ADD_MN_CASE_3X3( colwise_reduce, logsum, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_logsum )
	ADD_MN_CASE_3X3( rowwise_reduce, logsum, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_entropy )
	ADD_MN_CASE_3X3( colwise_reduce, entropy, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_entropy )
	ADD_MN_CASE_3X3( rowwise_reduce, entropy, DM, DN )
END_TPACK


BEGIN_TPACK( colwise_dot )
	ADD_MN_CASE_3X3( colwise_reduce, dot, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_dot )
	ADD_MN_CASE_3X3( rowwise_reduce, dot, DM, DN )
END_TPACK

BEGIN_TPACK( colwise_nrmdot )
	ADD_MN_CASE_3X3( colwise_reduce, nrmdot, DM, DN )
END_TPACK

BEGIN_TPACK( rowwise_nrmdot )
	ADD_MN_CASE_3X3( rowwise_reduce, nrmdot, DM, DN )
END_TPACK


BEGIN_MAIN_SUITE

	ADD_TPACK( colwise_sum )
	ADD_TPACK( rowwise_sum )

	ADD_TPACK( colwise_mean )
	ADD_TPACK( rowwise_mean )

	ADD_TPACK( colwise_maximum )
	ADD_TPACK( rowwise_maximum )

	ADD_TPACK( colwise_minimum )
	ADD_TPACK( rowwise_minimum )


	ADD_TPACK( colwise_L1norm )
	ADD_TPACK( rowwise_L1norm )

	ADD_TPACK( colwise_sqL2norm )
	ADD_TPACK( rowwise_sqL2norm )

	ADD_TPACK( colwise_L2norm )
	ADD_TPACK( rowwise_L2norm )

	ADD_TPACK( colwise_Linfnorm )
	ADD_TPACK( rowwise_Linfnorm )

	ADD_TPACK( colwise_logsum )
	ADD_TPACK( rowwise_logsum )

	ADD_TPACK( colwise_entropy )
	ADD_TPACK( rowwise_entropy )

	ADD_TPACK( colwise_dot )
	ADD_TPACK( rowwise_dot )

	ADD_TPACK( colwise_nrmdot )
	ADD_TPACK( rowwise_nrmdot )

END_MAIN_SUITE


