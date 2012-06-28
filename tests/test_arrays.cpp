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


typedef class lmat::darray<int, monitored_allocator<int> > blk_t;
typedef class lmat::scoped_array<int, monitored_allocator<int> > scblk_t;

lmat::test::memory_allocation_monitor lmat::test::global_memory_allocation_monitor;

#define ASSERT_PENDING(k) ASSERT_EQ(k, global_memory_allocation_monitor.num_pending_sections() )

#define ASSERT_NO_PENDING ASSERT_FALSE( global_memory_allocation_monitor.has_pending() );

SIMPLE_CASE( darray, constructs )
{
	ASSERT_NO_PENDING

	{
		blk_t a0;
		ASSERT_PENDING( 0 );

		ASSERT_EQ( a0.nelems(), 0 );
		ASSERT_EQ( a0.ptr_begin(), 0 );
		ASSERT_EQ( a0.ptr_end(), 0 );

		const index_t n = 5;
		blk_t a1(n);
		ASSERT_PENDING( 1 );

		ASSERT_EQ( a1.nelems(), n );
		ASSERT_NE( a1.ptr_begin(), 0 );
		ASSERT_EQ( a1.ptr_end(), a1.ptr_begin() + n );

		const int src[n] = {3, 4, 5, 6, 7};
		blk_t a2(n, src);
		ASSERT_PENDING( 2 );

		ASSERT_EQ( a2.nelems(), n );
		ASSERT_NE( a2.ptr_begin(), 0 );
		ASSERT_NE( a2.ptr_begin(), src );
		ASSERT_EQ( a2.ptr_end(), a2.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a2, src );

		int v = 8;
		const int csrc[n] = {v, v, v, v, v};
		blk_t a3(n, v);
		ASSERT_PENDING( 3 );

		ASSERT_EQ( a3.nelems(), n );
		ASSERT_NE( a3.ptr_begin(), 0 );
		ASSERT_EQ( a3.ptr_end(), a3.ptr_begin() + n );

		ASSERT_VEC_EQ( n, a3, csrc );
	}

	ASSERT_NO_PENDING
}


BEGIN_TPACK( darray )
	ADD_SIMPLE_CASE( darray, constructs )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( darray )
END_MAIN_SUITE




