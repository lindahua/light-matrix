/**
 * @file test_array.cpp
 *
 * Unit testing for array classes.
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/common/block.h>
#include "mon_alloc.h"

using namespace lmat;
using namespace lmat::test;

// explicit instantiation

template class lmat::dblock<double, aligned_allocator<double> >;
template class lmat::sblock<double, 4>;

typedef class lmat::dblock<int, monitored_allocator<int> > arr_t;

lmat::test::memory_allocation_monitor lmat::test::global_memory_allocation_monitor;

#define ASSERT_PENDING(k) ASSERT_EQ(k, global_memory_allocation_monitor.num_pending_sections() )

#define ASSERT_NO_PENDING ASSERT_FALSE( global_memory_allocation_monitor.has_pending() );


SIMPLE_CASE( dblock, construct )
{
	ASSERT_NO_PENDING

	{
		arr_t a0;
		ASSERT_PENDING( 0 );

		ASSERT_EQ( a0.nelems(), 0 );
		ASSERT_EQ( a0.ptr_data(), 0 );

		const index_t n = 5;
		arr_t a1(n);
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_data(), 0 );

		const int src[n] = {3, 4, 5, 6, 7};
		arr_t a2(n, copy_from(src));
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_data(), 0 );
		ASSERT_NE( a2.ptr_data(), src );

		ASSERT_VEC_EQ( n, a2, src );

		int v = 8;
		const int csrc[n] = {v, v, v, v, v};
		arr_t a3(n, fill(v));
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_data(), 0 );

		ASSERT_VEC_EQ( n, a3, csrc );

		const int zsrc[n] = {0, 0, 0, 0, 0};
		arr_t a4(n, zero<int>());
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a4.nelems(), n );
		ASSERT_NE( a4.ptr_data(), 0 );

		ASSERT_VEC_EQ( n, a4, zsrc );
	}

	ASSERT_NO_PENDING
}

SIMPLE_CASE( dblock, copy_and_assign )
{
	ASSERT_NO_PENDING

	{
		const index_t n = 5;
		const index_t n2 = 3;
		const int src[n] = {3, 4, 5, 6, 7};

		arr_t a1(n, copy_from(src));
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_data(), 0 );

		ASSERT_VEC_EQ( n, a1, src );

		arr_t a2(a1);
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_data(), 0 );
		ASSERT_NE( a2.ptr_data(), a1.ptr_data() );

		ASSERT_VEC_EQ( n, a2, src );

		arr_t a3;
		ASSERT_PENDING( 2 );

		a3 = a1;
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_data(), 0 );
		ASSERT_NE( a3.ptr_data(), a1.ptr_data() );

		ASSERT_VEC_EQ( n, a3, src );

		arr_t a4(n);
		ASSERT_PENDING( 4 );
		const int *p4_0 = a4.ptr_data();

		a4 = a1;
		ASSERT_PENDING( 4 );

		ASSERT_EQ( a4.nelems(), n );
		ASSERT_EQ( a4.ptr_data(), p4_0 );

		ASSERT_VEC_EQ( n, a4, src );

		arr_t a5(n2);
		ASSERT_PENDING( 5 );

		a5 = a1;
		ASSERT_PENDING( 5 );

		ASSERT_EQ( a5.nelems(), n );
		ASSERT_NE( a5.ptr_data(), 0 );
		ASSERT_NE( a5.ptr_data(), a1.ptr_data() );

		ASSERT_VEC_EQ( n, a5, src );
	}

	ASSERT_NO_PENDING
}


SIMPLE_CASE( dblock, swap )
{
	ASSERT_NO_PENDING

	{
		const index_t n1 = 5;
		const index_t n2 = 3;

		const int src1[n1] = {3, 4, 5, 6, 7};
		const int src2[n2] = {9, 8, 7};

		arr_t a0;
		arr_t a1(n1, copy_from(src1));
		const int *pa1 = a1.ptr_data();

		ASSERT_PENDING(1);

		ASSERT_EQ( a0.nelems(), 0 );
		ASSERT_EQ( a1.nelems(), n1 );
		ASSERT_VEC_EQ( n1, a1, src1 );

		swap(a0, a1);

		ASSERT_PENDING(1);

		ASSERT_EQ( a0.nelems(), n1 );
		ASSERT_EQ( a0.ptr_data(), pa1 );
		ASSERT_EQ( a1.nelems(), 0 );
		ASSERT_EQ( a1.ptr_data(), 0 );

		arr_t a2(n2, copy_from(src2));
		const int *pa2 = a2.ptr_data();

		ASSERT_PENDING(2);

		ASSERT_EQ( a2.nelems(), n2 );
		ASSERT_VEC_EQ( n2, a2, src2 );

		swap(a0, a2);

		ASSERT_PENDING(2);

		ASSERT_EQ( a0.nelems(), n2 );
		ASSERT_EQ( a0.ptr_data(), pa2 );
		ASSERT_EQ( a2.nelems(), n1 );
		ASSERT_EQ( a2.ptr_data(), pa1 );
	}

	ASSERT_NO_PENDING
}


SIMPLE_CASE( sblock, construct )
{
	const index_t n = 5;
	typedef sblock<int, n> sarr_t;

	sarr_t a1;

	ASSERT_EQ( a1.nelems(), n );
	ASSERT_NE( a1.ptr_data(), 0 );

	const int src[n] = {3, 4, 5, 6, 7};
	sarr_t a2(copy_from(src));

	ASSERT_EQ( a2.nelems(), n );
	ASSERT_NE( a2.ptr_data(), 0 );
	ASSERT_NE( a2.ptr_data(), src );

	ASSERT_VEC_EQ( n, a2, src );

	int v = 8;
	const int csrc[n] = {v, v, v, v, v};
	sarr_t a3(fill(v));

	ASSERT_EQ( a3.nelems(), n );
	ASSERT_NE( a3.ptr_data(), 0 );

	ASSERT_VEC_EQ( n, a3, csrc );
}


SIMPLE_CASE( sblock, copy_and_assign )
{
	const index_t n = 5;
	typedef sblock<int, n> sarr_t;

	const int src[n] = {3, 4, 5, 6, 7};

	sarr_t a1(copy_from(src));

	ASSERT_EQ( a1.nelems(), n );
	ASSERT_NE( a1.ptr_data(), 0 );

	ASSERT_VEC_EQ( n, a1, src );

	sarr_t a2(a1);

	ASSERT_EQ( a2.nelems(), n );
	ASSERT_NE( a2.ptr_data(), 0 );
	ASSERT_NE( a2.ptr_data(), a1.ptr_data() );

	ASSERT_VEC_EQ( n, a2, src );

	sarr_t a3;
	a3 = a1;

	ASSERT_EQ( a3.nelems(), n );
	ASSERT_NE( a3.ptr_data(), 0 );
	ASSERT_NE( a3.ptr_data(), a1.ptr_data() );

	ASSERT_VEC_EQ( n, a3, src );
}


SIMPLE_CASE( sblock, swap )
{
	const index_t n = 5;
	typedef sblock<int, n> sarr_t;

	const int src1[n] = {3, 4, 5, 6, 7};
	const int src2[n] = {9, 8, 7, 6, 5};

	sarr_t a1(copy_from(src1));
	sarr_t a2(copy_from(src2));

	ASSERT_VEC_EQ( n, a1, src1 );
	ASSERT_VEC_EQ( n, a2, src2 );

	swap(a1, a2);

	ASSERT_VEC_EQ( n, a1, src2 );
	ASSERT_VEC_EQ( n, a2, src1 );
}


BEGIN_TPACK( dblock )
	ADD_SIMPLE_CASE( dblock, construct )
	ADD_SIMPLE_CASE( dblock, copy_and_assign )
	ADD_SIMPLE_CASE( dblock, swap )
END_TPACK

BEGIN_TPACK( sblock )
	ADD_SIMPLE_CASE( sblock, construct )
	ADD_SIMPLE_CASE( sblock, copy_and_assign )
	ADD_SIMPLE_CASE( sblock, swap )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( dblock )
	ADD_TPACK( sblock )
END_MAIN_SUITE




