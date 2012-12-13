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

	static void set(simd_pack<float, sse_t>& pk, const float *s)
	{
		pk.set(s[0], s[1], s[2], s[3]);
	}
};

template<> struct elemwise_construct<double>
{
	static simd_pack<double, sse_t> get(const double* s)
	{
		return simd_pack<double, sse_t>(s[0], s[1]);
	}

	static void set(simd_pack<double, sse_t>& pk, const double *s)
	{
		pk.set(s[0], s[1]);
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


T_CASE( sse_pack, sets )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	pack_t pk;

	T v1 = T(3.2);
	pk.set(v1);
	ASSERT_SIMD_EQ( pk, v1 );

	T r2[width];
	for (unsigned i = 0; i < width; ++i) r2[i] = T(2.5 + i);
	elemwise_construct<T>::set(pk, r2);
	ASSERT_SIMD_EQ(pk, r2);

	T v0 = T(0);
	pk.reset();
	ASSERT_SIMD_EQ(pk, v0);
}


T_CASE( sse_pack, loads )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	const unsigned int len = 2 * width + 1;
	LMAT_ALIGN_SSE T src[len];
	for (unsigned i = 0; i < len; ++i) src[i] = T(1.8 + i);

	pack_t pk = pack_t::zeros();

	pk.load_a(src);
	ASSERT_SIMD_EQ(pk, src);

	pk.load_u(src + 1);
	ASSERT_SIMD_EQ(pk, src + 1);
}

T_CASE( sse_pack, stores )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_SSE T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(1.8 + i);

	const unsigned int len = 2 * width + 1;
	LMAT_ALIGN_SSE T dst[len];

	pack_t pk;
	pk.load_a(src);

	for (unsigned i = 0; i < len; ++i) dst[i] = T(0);
	pk.store_a(dst);
	ASSERT_VEC_EQ(width, dst, src);

	for (unsigned i = 0; i < len; ++i) dst[i] = T(0);
	pk.store_u(dst + 1);
	ASSERT_VEC_EQ(width, dst + 1, src);
}


TI_CASE( sse_pack, load_parts )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_SSE T src_base[width + 1];
	T *src = src_base + 1;
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_part(I, src);

	T r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = T(0);
	for (int i = 0; i < I; ++i) r[i] = src[i];

	ASSERT_SIMD_EQ( pk, r );
}

TI_CASE( sse_pack, store_parts )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_SSE T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	T v = T(2.3);
	T r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = v;
	for (int i = 0; i < I; ++i) r[i] = src[i];

	LMAT_ALIGN_SSE T dst_base[width + 1];
	T *dst = dst_base + 1;
	for (unsigned i = 0; i < width; ++i) dst[i] = v;

	pk.store_part(I, dst);
	ASSERT_VEC_EQ( width, dst, r );
}

T_CASE( sse_pack, to_scalar )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_SSE T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	T v = pk.to_scalar();

	ASSERT_EQ(v, src[0]);
}

T_CASE( sse_pack, extracts )
{
	typedef simd_pack<T, sse_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_SSE T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	for (int i = 0; i < (int)width; ++i)
	{
		T v = pk.extract(i);
		ASSERT_EQ(v, src[i]);
	}
}

DEF_TPACK( sse_pack, constructs )
DEF_TPACK( sse_pack, sets )
DEF_TPACK( sse_pack, loads )
DEF_TPACK( sse_pack, stores )

BEGIN_TPACK( sse_pack_load_parts )
	ADD_TI_CASE( sse_pack, load_parts, float, 0 )
	ADD_TI_CASE( sse_pack, load_parts, float, 1 )
	ADD_TI_CASE( sse_pack, load_parts, float, 2 )
	ADD_TI_CASE( sse_pack, load_parts, float, 3 )
	ADD_TI_CASE( sse_pack, load_parts, float, 4 )

	ADD_TI_CASE( sse_pack, load_parts, double, 0 )
	ADD_TI_CASE( sse_pack, load_parts, double, 1 )
	ADD_TI_CASE( sse_pack, load_parts, double, 2 )
END_TPACK

BEGIN_TPACK( sse_pack_store_parts )
	ADD_TI_CASE( sse_pack, store_parts, float, 0 )
	ADD_TI_CASE( sse_pack, store_parts, float, 1 )
	ADD_TI_CASE( sse_pack, store_parts, float, 2 )
	ADD_TI_CASE( sse_pack, store_parts, float, 3 )
	ADD_TI_CASE( sse_pack, store_parts, float, 4 )

	ADD_TI_CASE( sse_pack, store_parts, double, 0 )
	ADD_TI_CASE( sse_pack, store_parts, double, 1 )
	ADD_TI_CASE( sse_pack, store_parts, double, 2 )
END_TPACK

DEF_TPACK( sse_pack, to_scalar )
DEF_TPACK( sse_pack, extracts )


BEGIN_MAIN_SUITE
	ADD_TPACK( sse_pack_constructs )
	ADD_TPACK( sse_pack_sets )
	ADD_TPACK( sse_pack_loads )
	ADD_TPACK( sse_pack_stores )
	ADD_TPACK( sse_pack_load_parts )
	ADD_TPACK( sse_pack_store_parts )
	ADD_TPACK( sse_pack_to_scalar )
	ADD_TPACK( sse_pack_extracts )
END_MAIN_SUITE

