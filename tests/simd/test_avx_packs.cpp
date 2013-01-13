/**
 * @file test_avx_packs.cpp
 *
 * @brief Unit tests of AVX packs
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/simd/avx_packs.h>
#include <cmath>

using namespace lmat;
using namespace lmat::test;

static_assert(simd_pack<float,  avx_t>::pack_width == 8, "Unexpected pack width");
static_assert(simd_pack<double, avx_t>::pack_width == 4, "Unexpected pack width");

template<typename T> struct elemwise_construct;

template<> struct elemwise_construct<float>
{
	static simd_pack<float, avx_t> get(const float* s)
	{
		return simd_pack<float, avx_t>(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
	}

	static void set(simd_pack<float, avx_t>& pk, const float *s)
	{
		pk.set(s[0], s[1], s[2], s[3], s[4], s[5], s[6], s[7]);
	}
};

template<> struct elemwise_construct<double>
{
	static simd_pack<double, avx_t> get(const double* s)
	{
		return simd_pack<double, avx_t>(s[0], s[1], s[2], s[3]);
	}

	static void set(simd_pack<double, avx_t>& pk, const double *s)
	{
		pk.set(s[0], s[1], s[2], s[3]);
	}
};


T_CASE( avx_pack_constructs )
{
	typedef simd_pack<T, avx_t> pack_t;
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

	pack_t pk3(r2);
	ASSERT_SIMD_EQ( pk3, r2 );

	pack_t pv1 = pack_t::ones();
	ASSERT_SIMD_EQ( pv1, T(1) );

	pack_t pv_inf = pack_t::inf();
	for (unsigned i = 0; i < width; ++i)
	{
		bool is_inf_i = std::isinf(pv_inf[i]) && pv_inf[i] > T(0);
		ASSERT_TRUE( is_inf_i );
	}

	pack_t pv_neginf = pack_t::neg_inf();
	for (unsigned i = 0; i < width; ++i)
	{
		bool is_neginf_i = std::isinf(pv_neginf[i]) && pv_neginf[i] < T(0);
		ASSERT_TRUE( is_neginf_i );
	}

	pack_t pv_nan = pack_t::nan();
	for (unsigned i = 0; i < width; ++i)
	{
		bool is_nan_i = std::isnan(pv_nan[i]);
		ASSERT_TRUE( is_nan_i );
	}
}



T_CASE( avx_pack_sets )
{
	typedef simd_pack<T, avx_t> pack_t;
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


T_CASE( avx_pack_loads )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	const unsigned int len = 2 * width + 1;
	LMAT_ALIGN_AVX T src[len];
	for (unsigned i = 0; i < len; ++i) src[i] = T(1.8 + i);

	pack_t pk = pack_t::zeros();

	pk.load_a(src);
	ASSERT_SIMD_EQ(pk, src);

	pk.load_u(src + 1);
	ASSERT_SIMD_EQ(pk, src + 1);
}

T_CASE( avx_pack_stores )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_AVX T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(1.8 + i);

	const unsigned int len = 2 * width + 1;
	LMAT_ALIGN_AVX T dst[len];

	pack_t pk;
	pk.load_a(src);

	for (unsigned i = 0; i < len; ++i) dst[i] = T(0);
	pk.store_a(dst);
	ASSERT_VEC_EQ(width, dst, src);

	for (unsigned i = 0; i < len; ++i) dst[i] = T(0);
	pk.store_u(dst + 1);
	ASSERT_VEC_EQ(width, dst + 1, src);
}


TI_CASE( avx_pack_load_parts )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_AVX T src_base[width + 1];
	T *src = src_base + 1;
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_part(I, src);

	T r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = T(0);
	for (int i = 0; i < I; ++i) r[i] = src[i];

	ASSERT_SIMD_EQ( pk, r );
}

TI_CASE( avx_pack_store_parts )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_AVX T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	T v = T(2.3);
	T r[width];
	for (unsigned i = 0; i < width; ++i) r[i] = v;
	for (int i = 0; i < I; ++i) r[i] = src[i];

	LMAT_ALIGN_AVX T dst_base[width + 1];
	T *dst = dst_base + 1;
	for (unsigned i = 0; i < width; ++i) dst[i] = v;

	pk.store_part(I, dst);
	ASSERT_VEC_EQ( width, dst, r );
}

T_CASE( avx_pack_to_scalar )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_AVX T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	T v = pk.to_scalar();

	ASSERT_EQ(v, src[0]);
}

T_CASE( avx_pack_extracts )
{
	typedef simd_pack<T, avx_t> pack_t;
	const unsigned int width = pack_t::pack_width;

	LMAT_ALIGN_AVX T src[width];
	for (unsigned i = 0; i < width; ++i) src[i] = T(2.4 + i);

	pack_t pk;
	pk.load_a(src);

	for (unsigned int i = 0; i < width; ++i)
	{
		T v = pk.extract(i);
		ASSERT_EQ(v, src[i]);
	}
}


AUTO_TPACK( avx_basics )
{
	ADD_T_CASE_FP( avx_pack_constructs )
	ADD_T_CASE_FP( avx_pack_sets )
	ADD_T_CASE_FP( avx_pack_loads )
	ADD_T_CASE_FP( avx_pack_stores )
}

AUTO_TPACK( avx_parts )
{
	ADD_TI_CASE( avx_pack_load_parts, float, 0 )
	ADD_TI_CASE( avx_pack_load_parts, float, 1 )
	ADD_TI_CASE( avx_pack_load_parts, float, 2 )
	ADD_TI_CASE( avx_pack_load_parts, float, 3 )
	ADD_TI_CASE( avx_pack_load_parts, float, 4 )

	ADD_TI_CASE( avx_pack_load_parts, double, 0 )
	ADD_TI_CASE( avx_pack_load_parts, double, 1 )
	ADD_TI_CASE( avx_pack_load_parts, double, 2 )

	ADD_TI_CASE( avx_pack_store_parts, float, 0 )
	ADD_TI_CASE( avx_pack_store_parts, float, 1 )
	ADD_TI_CASE( avx_pack_store_parts, float, 2 )
	ADD_TI_CASE( avx_pack_store_parts, float, 3 )
	ADD_TI_CASE( avx_pack_store_parts, float, 4 )

	ADD_TI_CASE( avx_pack_store_parts, double, 0 )
	ADD_TI_CASE( avx_pack_store_parts, double, 1 )
	ADD_TI_CASE( avx_pack_store_parts, double, 2 )
}

AUTO_TPACK( avx_elems )
{
	ADD_T_CASE_FP( avx_pack_to_scalar )
	ADD_T_CASE_FP( avx_pack_extracts )
}

