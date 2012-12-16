/**
 * @file test_avx_bpacks.cpp
 *
 * Unit testing of AVX boolean packs
 * 
 * @author Dahua Lin 
 */


#include "simd_test_base.h"
#include <light_mat/math/avx_bpacks.h>

using namespace lmat;
using namespace lmat::test;

using lmat::math::simd_bpack;
using lmat::math::avx_t;

static_assert(simd_bpack<float,  avx_t>::pack_width == 8, "Unexpected pack width");
static_assert(simd_bpack<double, avx_t>::pack_width == 4, "Unexpected pack width");

template<typename T> struct elemwise_construct;

template<> struct elemwise_construct<float>
{
	static simd_bpack<float, avx_t> get(const bool* s)
	{
		return simd_bpack<float, avx_t>(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
	}

	static void set(simd_bpack<float, avx_t>& pk, const bool *s)
	{
		pk.set(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
	}
};

template<> struct elemwise_construct<double>
{
	static simd_bpack<double, avx_t> get(const bool* s)
	{
		return simd_bpack<double, avx_t>(s[0], s[1], s[2], s[3]);
	}

	static void set(simd_bpack<double, avx_t>& pk, const bool *s)
	{
		pk.set(s[0], s[1], s[2], s[3]);
	}
};


#define DEF_TPACK( pname, tname ) \
	BEGIN_TPACK( pname##_##tname ) \
		ADD_T_CASE( pname, tname, float ) \
		ADD_T_CASE( pname, tname, double ) \
	END_TPACK


T_CASE( avx_bpack, constructs )
{
	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;
	const unsigned int width = bpack_t::pack_width;

	bpack_t pk0 = bpack_t::all_false();
	ASSERT_SIMD_EQ( pk0,  bint(0));

	bpack_t pk1 = bpack_t::all_true();
	ASSERT_SIMD_EQ( pk1,  bint(-1));

	bpack_t pk2( false );
	ASSERT_SIMD_EQ( pk2, bint(0) );

	bpack_t pk3( true );
	ASSERT_SIMD_EQ( pk3, bint(-1) );

	bool s[width];
	for (unsigned i = 0; i < width; ++i) s[i] = (i % 2 == 0);

	bint r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = -bint(s[i]);

	bpack_t pk4 = elemwise_construct<T>::get(s);
	ASSERT_SIMD_EQ( pk4, r );
}


T_CASE( avx_bpack, load_and_store )
{
	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;
	const unsigned int width = bpack_t::pack_width;

	bool s[width];
	bint si[width];
	bool r[width];

	for (unsigned i = 0; i < width; ++i)
	{
		s[i] = (i % 2 == 0);
		si[i] = -bint(s[i]);
		r[i] = false;
	}

	bpack_t pk(s);
	ASSERT_SIMD_EQ( pk, si );

	pk.store(r);
	ASSERT_VEC_EQ( width, s, r);
}


T_CASE( avx_bpack, set )
{
	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;
	const unsigned int width = bpack_t::pack_width;

	bpack_t pk;

	pk.set( true );
	ASSERT_SIMD_EQ( pk,  bint(-1));

	pk.set( false );
	ASSERT_SIMD_EQ( pk,  bint(0));

	bool s[width];
	for (unsigned i = 0; i < width; ++i) s[i] = (i % 2 == 0);

	bint r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = -bint(s[i]);

	elemwise_construct<T>::set(pk, s);
	ASSERT_SIMD_EQ( pk, r );
}


T_CASE( avx_bpack, to_scalar )
{
	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;
	const unsigned int width = bpack_t::pack_width;

	bpack_t pk;
	pk.set( true );
	ASSERT_EQ( pk.to_scalar(), true );

	pk.set( false );
	ASSERT_EQ( pk.to_scalar(), false );

	bool s[width];
	for (unsigned i = 0; i < width; ++i) s[i] = (i % 2 == 0);

	elemwise_construct<T>::set(pk, s);
	ASSERT_EQ( pk.to_scalar(), true );
}


T_CASE( avx_bpack, extracts )
{
	typedef simd_bpack<T, avx_t> bpack_t;
	typedef typename bpack_t::bint_type bint;
	const unsigned int width = bpack_t::pack_width;

	bool s[width];
	for (unsigned i = 0; i < width; ++i) s[i] = (i % 2 == 0);

	bpack_t pk;
	elemwise_construct<T>::set(pk, s);

	for (unsigned int i = 0; i < width; ++i)
	{
		ASSERT_EQ( pk.extract(i), s[i] );
	}

	for (unsigned i = 0; i < width; ++i) s[i] = (i % 3 == 0);

	elemwise_construct<T>::set(pk, s);

	for (unsigned int i = 0; i < width; ++i)
	{
		ASSERT_EQ( pk.extract(i), s[i] );
	}
}


DEF_TPACK( avx_bpack, constructs )
DEF_TPACK( avx_bpack, load_and_store )
DEF_TPACK( avx_bpack, set )
DEF_TPACK( avx_bpack, to_scalar )
DEF_TPACK( avx_bpack, extracts )


BEGIN_MAIN_SUITE
	ADD_TPACK( avx_bpack_constructs )
	ADD_TPACK( avx_bpack_load_and_store )
	ADD_TPACK( avx_bpack_set )
	ADD_TPACK( avx_bpack_to_scalar )
	ADD_TPACK( avx_bpack_extracts )
END_MAIN_SUITE



