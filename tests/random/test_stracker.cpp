/**
 * @file test_stracker.cpp
 *
 * @brief Unit testing of stream tracker
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/random/stream_tracker.h>

using namespace lmat;
using namespace lmat::random;
using namespace lmat::test;

typedef stream_tracker<uint32_t> trk_t;

SIMPLE_CASE( stracker, basics )
{
	size_t len = 20;
	trk_t t(len);

	ASSERT_EQ( t.length(), len );
	ASSERT_EQ( t.offset(), len );
	ASSERT_EQ( t.remain(), 0 );
	ASSERT_EQ( t.is_end(), true );

	t.rewind();

	ASSERT_EQ( t.length(), len );
	ASSERT_EQ( t.offset(), 0 );
	ASSERT_EQ( t.remain(), len );
	ASSERT_EQ( t.is_end(), false );

	t.set_end();

	ASSERT_EQ( t.length(), len );
	ASSERT_EQ( t.offset(), len );
	ASSERT_EQ( t.remain(), 0 );
	ASSERT_EQ( t.is_end(), true );
}


template<typename BTag>
void test_stracker_tobound(unsigned int w, size_t len)
{
	trk_t t(len);
	for (size_t i = 0; i <= len; ++i)
	{
		t.set_offset(i);
		t.to_boundary(BTag());

		size_t o = (i / w) * w;
		if (o < i) o += w;

		size_t r = o < len ? len - o : 0;

		ASSERT_EQ( t.offset(), o );
		ASSERT_EQ( t.remain(), r );
		ASSERT_EQ( t.is_end(), o >= len );
	}
}

SIMPLE_CASE( stracker, to_boundary_dbl )
{
	test_stracker_tobound<bdtags::dbl>(2, 20);
}

SIMPLE_CASE( stracker, to_boundary_quad )
{
	test_stracker_tobound<bdtags::quad>(4, 20);
}

SIMPLE_CASE( stracker, to_boundary_oct )
{
	test_stracker_tobound<bdtags::oct>(8, 30);
}

SIMPLE_CASE( stracker, to_boundary_hex )
{
	test_stracker_tobound<bdtags::hex>(16, 50);
}


SIMPLE_CASE( stracker, forward_bytes )
{
	size_t len = 10;
	trk_t t(len);

	for (size_t i = 0; i <= 10; ++i)
	{
		size_t rb = (len - i) * 4;
		for (size_t n = 0; n <= rb; ++n)
		{
			t.set_offset(i);
			t.forward_bytes(n);

			size_t o = i + n / 4;
			if (n % 4) ++o;

			ASSERT_EQ( t.offset(), o );
		}
	}
}


BEGIN_TPACK( stracker )
	ADD_SIMPLE_CASE( stracker, basics )
	ADD_SIMPLE_CASE( stracker, to_boundary_dbl )
	ADD_SIMPLE_CASE( stracker, to_boundary_quad )
	ADD_SIMPLE_CASE( stracker, to_boundary_oct )
	ADD_SIMPLE_CASE( stracker, to_boundary_hex )
	ADD_SIMPLE_CASE( stracker, forward_bytes )
END_TPACK


BEGIN_MAIN_SUITE
	ADD_TPACK( stracker )
END_MAIN_SUITE




