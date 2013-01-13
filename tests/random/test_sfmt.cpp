/**
 * @file test_sfmt.cpp
 *
 * @brief Test of SFMT stream
 *
 * @author Dahua Lin
 */

#include "../test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/random/sfmt.h>
#include <cstring>
#include <fstream>

using namespace lmat;
using namespace lmat::random;
using namespace lmat::test;


const index_t vlen = 1000;
const uint32_t seed0 = 1234;

dense_col<uint32_t> read_sfmt_vdata(unsigned int mexp)
{
	dense_col<uint32_t> data(vlen);

	char filename[128];
	std::snprintf(filename, 128, "data/sfmt/sfmt.%06u.txt", mexp);

	std::ifstream fin;
	fin.open(filename);

	if (fin.fail())
		throw std::runtime_error("SFMT data file open failed");

	for (index_t i = 0; i < vlen; ++i)
	{
		fin >> data[i];
	}

	fin.close();
/*
	for (index_t i = 0; i < vlen; ++i)
	{
		std::printf("%11u", data[i]);
		if ((i+1) % 5 == 0) std::printf("\n");
	}
*/
	return data;
}


template<unsigned int MEXP>
void verify_sfmt_stream()
{
	dense_col<uint32_t> vdata = read_sfmt_vdata(MEXP);

	sfmt_rand_stream<MEXP> rs;

	dense_col<uint32_t> vgen(vlen);
	for (index_t i = 0; i < vlen; ++i)
	{
		vgen[i] = rs.rand_u32();
	}
	ASSERT_VEC_EQ( vlen, vgen, vdata );
}

template<unsigned int MEXP>
void verify_sfmt_u64()
{
	sfmt_rand_stream<MEXP> rs;

	dense_col<uint32_t> v32(vlen);
	for (index_t i = 0; i < vlen; ++i)
	{
		v32[i] = rs.rand_u32();
	}

	rs.set_seed(seed0);

	index_t n = vlen / 2;

	for (index_t i = 0; i < n; ++i)
	{
		uint64_t x = rs.rand_u64();

		uint32_t l = v32[2 * i];
		uint32_t h = v32[2 * i + 1];
		uint64_t x0 = (uint64_t)l | ((uint64_t)h << 32);

		ASSERT_EQ( x, x0 );
	}

	rs.set_seed(seed0);
	rs.rand_u32(); // ignore one unit

	for (index_t i = 0; i < n-1; ++i)
	{
		uint64_t x = rs.rand_u64();
		uint32_t l = v32[2 * i + 2];
		uint32_t h = v32[2 * i + 3];
		uint64_t x0 = (uint64_t)l | ((uint64_t)h << 32);

		ASSERT_EQ( x, x0 );
	}
}


template<unsigned int MEXP>
void verify_sfmt_m128()
{
	sfmt_rand_stream<MEXP> rs;

	dense_col<uint32_t> v32(vlen);
	for (index_t i = 0; i < vlen; ++i)
	{
		v32[i] = rs.rand_u32();
	}

	rs.set_seed(seed0);

	index_t n = vlen / 4;

	LMAT_ALIGN(16) uint32_t x[4];

	for (index_t i = 0; i < n; ++i)
	{
		__m128i p = rs.rand_pack(sse_t());
		_mm_store_si128(reinterpret_cast<__m128i*>(x), p);

		const uint32_t *x0 = &v32[i * 4];

		ASSERT_VEC_EQ( 4, x, x0 );
	}

	for (index_t o = 1; o < 4; ++o)
	{
		rs.set_seed(seed0);
		for (index_t j = 0; j < o; ++j) rs.rand_u32(); // ignore o units

		for (index_t i = 0; i < n-1; ++i)
		{
			__m128i p = rs.rand_pack(sse_t());
			_mm_store_si128(reinterpret_cast<__m128i*>(x), p);

			const uint32_t *x0 = &v32[(i+1) * 4];

			ASSERT_VEC_EQ(4, x, x0);
		}
	}
}


#ifdef LMAT_HAS_AVX

template<unsigned int MEXP>
void verify_sfmt_m256()
{
	sfmt_rand_stream<MEXP> rs;

	dense_col<uint32_t> v32(vlen);
	for (index_t i = 0; i < vlen; ++i)
	{
		v32[i] = rs.rand_u32();
	}

	rs.set_seed(seed0);

	index_t n = vlen / 8;

	LMAT_ALIGN(32) uint32_t x[8];

	for (index_t i = 0; i < n; ++i)
	{
		__m256i p = rs.rand_pack(avx_t());
		_mm256_store_si256(reinterpret_cast<__m256i*>(x), p);

		const uint32_t *x0 = &v32[i * 8];

		ASSERT_VEC_EQ( 8, x, x0 );
	}

	for (index_t o = 1; o < 8; ++o)
	{
		rs.set_seed(seed0);
		for (index_t j = 0; j < o; ++j) rs.rand_u32(); // ignore o units

		for (index_t i = 0; i < n-1; ++i)
		{
			__m256i p = rs.rand_pack(avx_t());
			_mm256_store_si256(reinterpret_cast<__m256i*>(x), p);

			const uint32_t *x0 = &v32[(i+1) * 8];

			ASSERT_VEC_EQ(8, x, x0);
		}
	}
}

#endif


template<unsigned int MEXP>
void test_sfmt_randseq(sfmt_rand_stream<MEXP>& rs, index_t ignore, index_t n)  // n units
{
	dense_row<uint32_t> r(n);
	dense_row<uint32_t> x(n, zero());
	dense_row<uint32_t> r2(n);
	dense_row<uint32_t> x2(n, zero());

	// generate one by one

	rs.set_seed(seed0);
	for (index_t i = 0; i < ignore; ++i) rs.rand_u32();
	for (index_t i = 0; i < n; ++i) r[i] = rs.rand_u32();
	for (index_t i = 0; i < n; ++i) r2[i] = rs.rand_u32();

	// generate in batch

	rs.set_seed(seed0);
	for (index_t i = 0; i < ignore; ++i) rs.rand_u32();
	rs.rand_seq((size_t)n * sizeof(uint32_t), x.ptr_data());
	rs.rand_seq((size_t)n * sizeof(uint32_t), x2.ptr_data());

	// compare

	ASSERT_VEC_EQ(n, x, r);
	ASSERT_VEC_EQ(n, x2, r2);
}


template<unsigned int MEXP>
void verify_sfmt_seq()
{
	sfmt_rand_stream<MEXP> rs;
	index_t nu = (index_t)rs.internal_nunits();

	const int nstarts = 4;
	index_t starts[nstarts] = {0, 3, nu / 4, nu / 2};

	const int nlens = 6;
	index_t lens[nlens] = {3, nu / 4, nu / 2, nu, (nu * 2 + nu / 2), (nu * 3 + nu / 4) };

	for (int i = 0; i < nstarts; ++i)
	{
		for (int j = 0; j < nlens; ++j)
		{
			index_t start = starts[i];
			index_t len = lens[j];

			test_sfmt_randseq(rs, start, len);
		}
	}
}



#define DEF_SFMT_TESTS( packname, tfunname ) \
		SIMPLE_CASE( packname##_1279 ) { tfunname<1279>(); } \
		SIMPLE_CASE( packname##_2281 ) { tfunname<2281>(); } \
		SIMPLE_CASE( packname##_4253 ) { tfunname<4253>(); } \
		SIMPLE_CASE( packname##_11213 ) { tfunname<11213>(); } \
		SIMPLE_CASE( packname##_19937 ) { tfunname<19937>(); } \
		SIMPLE_CASE( packname##_44497 ) { tfunname<44497>(); } \
		SIMPLE_CASE( packname##_86243 ) { tfunname<86243>(); } \
		SIMPLE_CASE( packname##_132049 ) { tfunname<132049>(); } \
		AUTO_TPACK( packname ) { \
			ADD_SIMPLE_CASE( packname##_1279 ) \
			ADD_SIMPLE_CASE( packname##_2281 ) \
			ADD_SIMPLE_CASE( packname##_4253 ) \
			ADD_SIMPLE_CASE( packname##_11213 ) \
			ADD_SIMPLE_CASE( packname##_19937 ) \
			ADD_SIMPLE_CASE( packname##_44497 ) \
			ADD_SIMPLE_CASE( packname##_86243 ) \
			ADD_SIMPLE_CASE( packname##_132049 ) \
		}

DEF_SFMT_TESTS( sfmt_verify, verify_sfmt_stream )
DEF_SFMT_TESTS( sfmt_verify_u64, verify_sfmt_u64 )
DEF_SFMT_TESTS( sfmt_verify_m128, verify_sfmt_m128 )

#ifdef LMAT_HAS_AVX
DEF_SFMT_TESTS( sfmt_verify_m256, verify_sfmt_m256 )
#endif

DEF_SFMT_TESTS( sfmt_verify_seq, verify_sfmt_seq )


