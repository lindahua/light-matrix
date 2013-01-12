/**
 * @file test_sse_round.cpp
 *
 * @brief Unit testing of SSE rounding functions
 *
 * @author Dahua Lin
 */


#include "simd_test_base.h"
#include <light_mat/simd/sse_arith.h>

using namespace lmat;
using namespace lmat::test;

using lmat::math::floor;
using lmat::math::ceil;
using lmat::math::trunc;
using lmat::math::round;

using lmat::internal::floor_sse2;
using lmat::internal::ceil_sse2;
using lmat::internal::trunc_sse2;
using lmat::internal::round_sse2;

static float sa1[4] = { -1.2f, -1.5f, -1.7f, -2.0f };
static float sa2[4] = { 2.2f, 2.5f, 2.7f, 3.0f };

static double da1[2] = { -1.2, -1.5 };
static double da2[2] = { -1.7, -2.0 };
static double da3[2] = { 2.2, 2.5 };
static double da4[2] = { 2.7, 3.0 };


SIMPLE_CASE( sse_floor_f32 )
{
	sse_f32pk a1; a1.load_u(sa1);
	sse_f32pk a2; a2.load_u(sa2);

	float r1[4] = { -2.0f, -2.0f, -2.0f, -2.0f };
	float r2[4] = { 2.0f, 2.0f, 2.0f, 3.0f };

	ASSERT_SIMD_EQ( floor_sse2(a1), r1 );
	ASSERT_SIMD_EQ( floor_sse2(a2), r2 );

	ASSERT_SIMD_EQ( floor(a1), r1 );
	ASSERT_SIMD_EQ( floor(a2), r2 );
}

SIMPLE_CASE( sse_floor_f64 )
{
	sse_f64pk a1; a1.load_u(da1);
	sse_f64pk a2; a2.load_u(da2);
	sse_f64pk a3; a3.load_u(da3);
	sse_f64pk a4; a4.load_u(da4);

	double r1[2] = { -2.0, -2.0 };
	double r2[2] = { -2.0, -2.0 };
	double r3[2] = { 2.0, 2.0 };
	double r4[2] = { 2.0, 3.0 };

	ASSERT_SIMD_EQ( floor_sse2(a1), r1 );
	ASSERT_SIMD_EQ( floor_sse2(a2), r2 );
	ASSERT_SIMD_EQ( floor_sse2(a3), r3 );
	ASSERT_SIMD_EQ( floor_sse2(a4), r4 );

	ASSERT_SIMD_EQ( floor(a1), r1 );
	ASSERT_SIMD_EQ( floor(a2), r2 );
	ASSERT_SIMD_EQ( floor(a3), r3 );
	ASSERT_SIMD_EQ( floor(a4), r4 );
}


SIMPLE_CASE( sse_ceil_f32 )
{
	sse_f32pk a1; a1.load_u(sa1);
	sse_f32pk a2; a2.load_u(sa2);

	float r1[4] = { -1.0f, -1.0f, -1.0f, -2.0f };
	float r2[4] = { 3.0f, 3.0f, 3.0f, 3.0f };

	ASSERT_SIMD_EQ( ceil_sse2(a1), r1 );
	ASSERT_SIMD_EQ( ceil_sse2(a2), r2 );

	ASSERT_SIMD_EQ( ceil(a1), r1 );
	ASSERT_SIMD_EQ( ceil(a2), r2 );
}

SIMPLE_CASE( sse_ceil_f64 )
{
	sse_f64pk a1; a1.load_u(da1);
	sse_f64pk a2; a2.load_u(da2);
	sse_f64pk a3; a3.load_u(da3);
	sse_f64pk a4; a4.load_u(da4);

	double r1[2] = { -1.0, -1.0 };
	double r2[2] = { -1.0, -2.0 };
	double r3[2] = { 3.0, 3.0 };
	double r4[2] = { 3.0, 3.0 };

	ASSERT_SIMD_EQ( ceil_sse2(a1), r1 );
	ASSERT_SIMD_EQ( ceil_sse2(a2), r2 );
	ASSERT_SIMD_EQ( ceil_sse2(a3), r3 );
	ASSERT_SIMD_EQ( ceil_sse2(a4), r4 );

	ASSERT_SIMD_EQ( ceil(a1), r1 );
	ASSERT_SIMD_EQ( ceil(a2), r2 );
	ASSERT_SIMD_EQ( ceil(a3), r3 );
	ASSERT_SIMD_EQ( ceil(a4), r4 );
}


SIMPLE_CASE( sse_trunc_f32 )
{
	sse_f32pk a1; a1.load_u(sa1);
	sse_f32pk a2; a2.load_u(sa2);

	float r1[4] = { -1.0f, -1.0f, -1.0f, -2.0f };
	float r2[4] = { 2.0f, 2.0f, 2.0f, 3.0f };

	ASSERT_SIMD_EQ( trunc_sse2(a1), r1 );
	ASSERT_SIMD_EQ( trunc_sse2(a2), r2 );

	ASSERT_SIMD_EQ( trunc(a1), r1 );
	ASSERT_SIMD_EQ( trunc(a2), r2 );
}

SIMPLE_CASE( sse_trunc_f64 )
{
	sse_f64pk a1; a1.load_u(da1);
	sse_f64pk a2; a2.load_u(da2);
	sse_f64pk a3; a3.load_u(da3);
	sse_f64pk a4; a4.load_u(da4);

	double r1[2] = { -1.0, -1.0 };
	double r2[2] = { -1.0, -2.0 };
	double r3[2] = { 2.0, 2.0 };
	double r4[2] = { 2.0, 3.0 };

	ASSERT_SIMD_EQ( trunc_sse2(a1), r1 );
	ASSERT_SIMD_EQ( trunc_sse2(a2), r2 );
	ASSERT_SIMD_EQ( trunc_sse2(a3), r3 );
	ASSERT_SIMD_EQ( trunc_sse2(a4), r4 );

	ASSERT_SIMD_EQ( trunc(a1), r1 );
	ASSERT_SIMD_EQ( trunc(a2), r2 );
	ASSERT_SIMD_EQ( trunc(a3), r3 );
	ASSERT_SIMD_EQ( trunc(a4), r4 );
}


SIMPLE_CASE( sse_round_f32 )
{
	sse_f32pk a1; a1.load_u(sa1);
	sse_f32pk a2; a2.load_u(sa2);

	float r1[4] = { -1.0f, -2.0f, -2.0f, -2.0f };
	float r2[4] = { 2.0f, 2.0f, 3.0f, 3.0f };

	ASSERT_SIMD_EQ( round_sse2(a1), r1 );
	ASSERT_SIMD_EQ( round_sse2(a2), r2 );

	ASSERT_SIMD_EQ( round(a1), r1 );
	ASSERT_SIMD_EQ( round(a2), r2 );
}

SIMPLE_CASE( sse_round_f64 )
{
	sse_f64pk a1; a1.load_u(da1);
	sse_f64pk a2; a2.load_u(da2);
	sse_f64pk a3; a3.load_u(da3);
	sse_f64pk a4; a4.load_u(da4);

	double r1[2] = { -1.0, -2.0 };
	double r2[2] = { -2.0, -2.0 };
	double r3[2] = { 2.0, 2.0 };
	double r4[2] = { 3.0, 3.0 };

	ASSERT_SIMD_EQ( round_sse2(a1), r1 );
	ASSERT_SIMD_EQ( round_sse2(a2), r2 );
	ASSERT_SIMD_EQ( round_sse2(a3), r3 );
	ASSERT_SIMD_EQ( round_sse2(a4), r4 );

	ASSERT_SIMD_EQ( round(a1), r1 );
	ASSERT_SIMD_EQ( round(a2), r2 );
	ASSERT_SIMD_EQ( round(a3), r3 );
	ASSERT_SIMD_EQ( round(a4), r4 );
}


AUTO_TPACK( sse_round )
{
	ADD_SIMPLE_CASE( sse_floor_f32 )
	ADD_SIMPLE_CASE( sse_floor_f64 )
	ADD_SIMPLE_CASE( sse_ceil_f32 )
	ADD_SIMPLE_CASE( sse_ceil_f64 )
	ADD_SIMPLE_CASE( sse_trunc_f32 )
	ADD_SIMPLE_CASE( sse_trunc_f64 )
	ADD_SIMPLE_CASE( sse_round_f32 )
	ADD_SIMPLE_CASE( sse_round_f64 )
}

