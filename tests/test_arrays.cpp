/**
 * @file test_array.cpp
 *
 * Unit testing for array classes.
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/core/array.h>
#include "mon_alloc.h"

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::darray<double, aligned_allocator<double> >;
template class lmat::scoped_array<double, aligned_allocator<double> >;
template class lmat::sarray<double, 4>;


typedef class lmat::darray<int, monitored_allocator<int> > arr_t;
typedef class lmat::scoped_array<int, monitored_allocator<int> > scarr_t;

lmat::test::memory_allocation_monitor lmat::test::global_memory_allocation_monitor;

#define ASSERT_PENDING(k) ASSERT_EQ(k, global_memory_allocation_monitor.num_pending_sections() )

#define ASSERT_NO_PENDING ASSERT_FALSE( global_memory_allocation_monitor.has_pending() );

SIMPLE_CASE( darray, construct )
{
	ASSERT_NO_PENDING

	{
		arr_t a0;
		ASSERT_PENDING( 0 );

		ASSERT_EQ( a0.nelems(), 0 );
		ASSERT_EQ( a0.ptr_begin(), 0 );
		ASSERT_EQ( a0.ptr_end(), 0 );

		const index_t n = 5;
		arr_t a1(n);
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_begin(), 0 );
		ASSERT_EQ( a1.ptr_end(), a1.ptr_begin() + n );

		const int src[n] = {3, 4, 5, 6, 7};
		arr_t a2(n, src);
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_begin(), 0 );
		ASSERT_NE( a2.ptr_begin(), src );
		ASSERT_EQ( a2.ptr_end(), a2.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a2, src );

		int v = 8;
		const int csrc[n] = {v, v, v, v, v};
		arr_t a3(n, v);
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_begin(), 0 );
		ASSERT_EQ( a3.ptr_end(), a3.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a3, csrc );
	}

	ASSERT_NO_PENDING
}

SIMPLE_CASE( darray, copy_and_assign )
{
	ASSERT_NO_PENDING

	{
		const index_t n = 5;
		const index_t n2 = 3;
		const int src[n] = {3, 4, 5, 6, 7};

		arr_t a1(n, src);
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_begin(), 0 );

		ASSERT_VEC_EQ( n, a1, src );

		arr_t a2(a1);
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_begin(), 0 );
		ASSERT_NE( a2.ptr_begin(), a1.ptr_begin() );

		ASSERT_VEC_EQ( n, a2, src );

		arr_t a3;
		ASSERT_PENDING( 2 );

		a3 = a1;
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_begin(), 0 );
		ASSERT_NE( a3.ptr_begin(), a1.ptr_begin() );

		ASSERT_VEC_EQ( n, a3, src );

		arr_t a4(n);
		ASSERT_PENDING( 4 );
		const int *p4_0 = a4.ptr_begin();

		a4 = a1;
		ASSERT_PENDING( 4 );

		ASSERT_EQ( a4.nelems(), n );
		ASSERT_EQ( a4.ptr_begin(), p4_0 );

		ASSERT_VEC_EQ( n, a4, src );

		arr_t a5(n2);
		ASSERT_PENDING( 5 );

		a5 = a1;
		ASSERT_PENDING( 5 );

		ASSERT_EQ( a5.nelems(), n );
		ASSERT_NE( a5.ptr_begin(), 0 );
		ASSERT_NE( a5.ptr_begin(), a1.ptr_begin() );

		ASSERT_VEC_EQ( n, a5, src );
	}

	ASSERT_NO_PENDING
}


SIMPLE_CASE( darray, swap )
{
	ASSERT_NO_PENDING

	{
		const index_t n1 = 5;
		const index_t n2 = 3;

		const int src1[n1] = {3, 4, 5, 6, 7};
		const int src2[n2] = {9, 8, 7};

		arr_t a0;
		arr_t a1(n1, src1);
		const int *pa1 = a1.ptr_begin();

		ASSERT_PENDING(1);

		ASSERT_EQ( a0.nelems(), 0 );
		ASSERT_EQ( a1.nelems(), n1 );
		ASSERT_VEC_EQ( n1, a1, src1 );

		swap(a0, a1);

		ASSERT_PENDING(1);

		ASSERT_EQ( a0.nelems(), n1 );
		ASSERT_EQ( a0.ptr_begin(), pa1 );
		ASSERT_EQ( a1.nelems(), 0 );
		ASSERT_EQ( a1.ptr_begin(), 0 );

		arr_t a2(n2, src2);
		const int *pa2 = a2.ptr_begin();

		ASSERT_PENDING(2);

		ASSERT_EQ( a2.nelems(), n2 );
		ASSERT_VEC_EQ( n2, a2, src2 );

		swap(a0, a2);

		ASSERT_PENDING(2);

		ASSERT_EQ( a0.nelems(), n2 );
		ASSERT_EQ( a0.ptr_begin(), pa2 );
		ASSERT_EQ( a2.nelems(), n1 );
		ASSERT_EQ( a2.ptr_begin(), pa1 );
	}

	ASSERT_NO_PENDING
}


SIMPLE_CASE( scoped_array, construct )
{
	ASSERT_NO_PENDING

	{
		const index_t n = 5;
		scarr_t a1(n);
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_begin(), 0 );
		ASSERT_EQ( a1.ptr_end(), a1.ptr_begin() + n );

		const int src[n] = {3, 4, 5, 6, 7};
		scarr_t a2(n, src);
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_begin(), 0 );
		ASSERT_NE( a2.ptr_begin(), src );
		ASSERT_EQ( a2.ptr_end(), a2.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a2, src );

		int v = 8;
		const int csrc[n] = {v, v, v, v, v};
		scarr_t a3(n, v);
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_begin(), 0 );
		ASSERT_EQ( a3.ptr_end(), a3.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a3, csrc );
	}

	ASSERT_NO_PENDING
}


SIMPLE_CASE( sarray, construct )
{
	const index_t n = 5;
	typedef sarray<int, n> sarr_t;

	sarr_t a1;

	ASSERT_EQ( a1.nelems(), n );
	ASSERT_NE( a1.ptr_begin(), 0 );
	ASSERT_EQ( a1.ptr_end(), a1.ptr_begin() + n );

	const int src[n] = {3, 4, 5, 6, 7};
	sarr_t a2(src);

	ASSERT_EQ( a2.nelems(), n );
	ASSERT_NE( a2.ptr_begin(), 0 );
	ASSERT_NE( a2.ptr_begin(), src );
	ASSERT_EQ( a2.ptr_end(), a2.ptr_begin() + n );

	ASSERT_VEC_EQ( n, a2, src );

	int v = 8;
	const int csrc[n] = {v, v, v, v, v};
	sarr_t a3(v);

	ASSERT_EQ( a3.nelems(), n );
	ASSERT_NE( a3.ptr_begin(), 0 );
	ASSERT_EQ( a3.ptr_end(), a3.ptr_begin() + n );

	ASSERT_VEC_EQ( n, a3, csrc );
}


SIMPLE_CASE( sarray, copy_and_assign )
{
	const index_t n = 5;
	typedef sarray<int, n> sarr_t;

	const int src[n] = {3, 4, 5, 6, 7};

	sarr_t a1(src);

	ASSERT_EQ( a1.nelems(), n );
	ASSERT_NE( a1.ptr_begin(), 0 );

	ASSERT_VEC_EQ( n, a1, src );

	sarr_t a2(a1);

	ASSERT_EQ( a2.nelems(), n );
	ASSERT_NE( a2.ptr_begin(), 0 );
	ASSERT_NE( a2.ptr_begin(), a1.ptr_begin() );

	ASSERT_VEC_EQ( n, a2, src );

	sarr_t a3;
	a3 = a1;

	ASSERT_EQ( a3.nelems(), n );
	ASSERT_NE( a3.ptr_begin(), 0 );
	ASSERT_NE( a3.ptr_begin(), a1.ptr_begin() );

	ASSERT_VEC_EQ( n, a3, src );
}


SIMPLE_CASE( sarray, swap )
{
	const index_t n = 5;
	typedef sarray<int, n> sarr_t;

	const int src1[n] = {3, 4, 5, 6, 7};
	const int src2[n] = {9, 8, 7, 6, 5};

	sarr_t a1(src1);
	sarr_t a2(src2);

	ASSERT_VEC_EQ( n, a1, src1 );
	ASSERT_VEC_EQ( n, a2, src2 );

	swap(a1, a2);

	ASSERT_VEC_EQ( n, a1, src2 );
	ASSERT_VEC_EQ( n, a2, src1 );
}


SIMPLE_CASE( array_op, copy )
{
	const index_t n = 5;
	const int src[n] = {3, 4, 5, 6, 7};
	int dst[n] = {0, 0, 0, 0, 0};

	arr_t a(n, -1);

	copy(src, a);

	ASSERT_VEC_EQ( n, a, src );

	copy(a, dst);

	ASSERT_VEC_EQ( n, dst, src );
}

SIMPLE_CASE( array_op, fill )
{
	const int n = 5;
	const int v = 3;
	const int r[n] = {v, v, v, v, v};

	arr_t a(n);
	fill(a, v);

	ASSERT_VEC_EQ( n, a, r );
}


SIMPLE_CASE( array_op, zero )
{
	const int n = 5;
	const int r[n] = {0, 0, 0, 0, 0};

	arr_t a(n, -1);
	zero(a);

	ASSERT_VEC_EQ( n, a, r );
}

SIMPLE_CASE( array_op, compare )
{
	const int n = 5;
	const int src1[n] = {1, 2, 3, 4, 5};
	const int src2[n] = {1, 2, 3, 4, 6};

	arr_t a1(n, src1);
	arr_t a2(n, src1);
	arr_t a3(n-1, src1);
	arr_t a4(n, src2);

	ASSERT_TRUE( a1 == a2 );
	ASSERT_FALSE( a1 != a2 );

	ASSERT_TRUE( a1 != a3 );
	ASSERT_FALSE( a1 == a3 );

	ASSERT_TRUE( a1 != a4 );
	ASSERT_FALSE( a1 == a4 );

	ASSERT_TRUE( elems_equal(a1, src1) );
	ASSERT_FALSE( elems_equal(a1, src2) );

	const int v = 7;
	ASSERT_FALSE( elems_equal(a2, v) );
	fill(a2, v);
	ASSERT_TRUE( elems_equal(a2, v) );
}



BEGIN_TPACK( darray )
	ADD_SIMPLE_CASE( darray, construct )
	ADD_SIMPLE_CASE( darray, copy_and_assign )
	ADD_SIMPLE_CASE( darray, swap )
END_TPACK


BEGIN_TPACK( scoped_array )
	ADD_SIMPLE_CASE( scoped_array, construct )
END_TPACK

BEGIN_TPACK( sarray )
	ADD_SIMPLE_CASE( sarray, construct )
	ADD_SIMPLE_CASE( sarray, copy_and_assign )
	ADD_SIMPLE_CASE( sarray, swap )
END_TPACK

BEGIN_TPACK( array_op )
	ADD_SIMPLE_CASE( array_op, copy )
	ADD_SIMPLE_CASE( array_op, fill )
	ADD_SIMPLE_CASE( array_op, zero )
	ADD_SIMPLE_CASE( array_op, compare )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( darray )
	ADD_TPACK( scoped_array )
	ADD_TPACK( sarray )
END_MAIN_SUITE




