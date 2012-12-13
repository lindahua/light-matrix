/**
 * @file test_sse_packs.cpp
 *
 * @brief Unit testing of SSE packs
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/math/sse_packs.h>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_pack;
using lmat::math::sse_t;

static_assert(simd_pack<float,  sse_t>::pack_width == 4, "Unexpected pack width");
static_assert(simd_pack<double, sse_t>::pack_width == 2, "Unexpected pack width");


template<typename T> struct elemwise_construct;

template<> struct elemwise_construct<float>
{
	static simd_pack<float, sse_t> get(const float* s)
	{
		return simd_pack<float, sse_t>(s[0], s[1], s[2], s[3]);
	}
};

template<> struct elemwise_construct<double>
{
	static simd_pack<double, sse_t> get(const double* s)
	{
		return simd_pack<double, sse_t>(s[0], s[1]);
	}
};


#define DEF_TPACK( pname, tname ) \
	BEGIN_TPACK( pname##_##tname ) \
		ADD_T_CASE( pname, tname, float ) \
		ADD_T_CASE( pname, tname, double ) \
	END_TPACK


T_CASE( sse_pack, constructs )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	pack_t pk0 = pack_t::zeros();
	ASSERT_EQ( pk0.width(), width );
	T v0 = T(0);
	ASSERT_SIMD_EQ( pk0, v0 );

	T v1 = T(2.5);
	pack_t pk1( v1 );
	ASSERT_SIMD_EQ( pk1, v1 );

	T r2[width];
	for (unsigned i = 0; i < width; ++i) r2[i] = T(1.5 + i);

	pack_t pk2 = elemwise_construct<T>::get(r2);
	ASSERT_SIMD_EQ( pk2, r2 );
}


DEF_TPACK( sse_pack, constructs )


BEGIN_MAIN_SUITE
	ADD_TPACK( sse_pack_constructs )
END_MAIN_SUITE

