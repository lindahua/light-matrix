/**
 * @file test_simd_vec.cpp
 *
 * @brief Unit testing of SIMD short vectors
 *
 * @author Dahua Lin
 */

#include "simd_test_base.h"
#include <light_mat/simd/simd_vec.h>
#include <cstring>

using namespace lmat;
using namespace lmat::test;


// Auxiliary functions

template<typename T, typename Kind, index_t Len>
struct svec_set_test;

template<typename T, typename Kind>
struct svec_set_test<T, Kind, 2>
{
	static void run()
	{
		typedef simd_vec<T, Kind, 2> vec_t;
		vec_t vec;
		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		const index_t len = 2;

		ASSERT_EQ( vec.length(), len );
		ASSERT_TRUE( inlen >= len );

		T vs[] = { T(1), T(2) };
		const T* pvec = (const T*)(&vec);

		// set same

		T r1[inlen];
		for (index_t i = 0; i < inlen; ++i) r1[i] = T(0);
		for (index_t i = 0; i < len; ++i) r1[i] = vs[0];

		vec.set(vs[0]);
		ASSERT_VEC_EQ( inlen, pvec, r1 );

		// set different

		T r2[inlen];
		for (index_t i = 0; i < inlen; ++i) r2[i] = T(0);
		for (index_t i = 0; i < len; ++i) r2[i] = vs[i];

		vec.set(vs[0], vs[1]);
		ASSERT_VEC_EQ( inlen, pvec, r2 );

		// reset to zero

		vec.reset();

		T r0[inlen];
		for (index_t i = 0; i < inlen; ++i) r0[i] = T(0);

		ASSERT_VEC_EQ( inlen, pvec, r0 );
	}
};

template<typename T, typename Kind>
struct svec_set_test<T, Kind, 3>
{
	static void run()
	{
		typedef simd_vec<T, Kind, 3> vec_t;
		vec_t vec;
		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		const index_t len = 3;

		ASSERT_EQ( vec.length(), len );
		ASSERT_TRUE( inlen >= len );

		T vs[] = { T(1), T(2), T(3) };
		const T* pvec = (const T*)(&vec);

		// set same

		T r1[inlen];
		for (index_t i = 0; i < inlen; ++i) r1[i] = T(0);
		for (index_t i = 0; i < len; ++i) r1[i] = vs[0];

		vec.set(vs[0]);
		ASSERT_VEC_EQ( inlen, pvec, r1 );

		// set different

		T r2[inlen];
		for (index_t i = 0; i < inlen; ++i) r2[i] = T(0);
		for (index_t i = 0; i < len; ++i) r2[i] = vs[i];

		vec.set(vs[0], vs[1], vs[2]);
		ASSERT_VEC_EQ( inlen, pvec, r2 );

		// reset to zero

		vec.reset();

		T r0[inlen];
		for (index_t i = 0; i < inlen; ++i) r0[i] = T(0);

		ASSERT_VEC_EQ( inlen, pvec, r0 );
	}
};

template<typename T, typename Kind>
struct svec_set_test<T, Kind, 4>
{
	static void run()
	{
		typedef simd_vec<T, Kind, 4> vec_t;
		vec_t vec;
		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		const index_t len = 4;

		ASSERT_EQ( vec.length(), len );
		ASSERT_TRUE( inlen >= len );

		T vs[] = { T(1), T(2), T(3), T(4) };
		const T* pvec = (const T*)(&vec);

		// set same

		T r1[inlen];
		for (index_t i = 0; i < inlen; ++i) r1[i] = T(0);
		for (index_t i = 0; i < len; ++i) r1[i] = vs[0];

		vec.set(vs[0]);
		ASSERT_VEC_EQ( inlen, pvec, r1 );

		// set different

		T r2[inlen];
		for (index_t i = 0; i < inlen; ++i) r2[i] = T(0);
		for (index_t i = 0; i < len; ++i) r2[i] = vs[i];

		vec.set(vs[0], vs[1], vs[2], vs[3]);
		ASSERT_VEC_EQ( inlen, pvec, r2 );

		// reset to zero

		vec.reset();

		T r0[inlen];
		for (index_t i = 0; i < inlen; ++i) r0[i] = T(0);

		ASSERT_VEC_EQ( inlen, pvec, r0 );
	}
};


TN_CASE( sse_svec_set )
{
	svec_set_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_set )
{
	svec_set_test<T, avx_t, N>::run();
}
#endif

template<typename T, typename Kind, index_t Len>
struct svec_load_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;
		vec_t vec;
		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));

		LMAT_ALIGN(32) T vs[Len + 2];
		for (index_t i = 0; i < Len + 2; ++i) vs[i] = T(2 * i + 1);
		const T* pvec = (const T*)(&vec);

		T r1[inlen];
		for (index_t i = 0; i < inlen; ++i) r1[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r1[i] = vs[i];

		vec.load_a(vs);
		ASSERT_VEC_EQ( inlen, pvec, r1 );

		T r2[inlen];
		for (index_t i = 0; i < inlen; ++i) r2[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r2[i] = vs[i + 1];

		vec.load_u(vs + 1);
		ASSERT_VEC_EQ( inlen, pvec, r2 );
	}
};


TN_CASE( sse_svec_load )
{
	svec_load_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_load )
{
	svec_load_test<T, avx_t, N>::run();
}
#endif

template<typename T, typename Kind, index_t Len>
struct svec_store_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs[Len];
		for (index_t i = 0; i < Len; ++i) vs[i] = T(2 * i + 1);
		vec_t vec(vs);

		LMAT_ALIGN(32) T dst[2 * Len];
		for (index_t i = 0; i < 2 * Len; ++i) dst[i] = T(0);
		vec.store_a(dst);

		T r[2 * Len];
		for (index_t i = 0; i < 2 * Len; ++i) r[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r[i] = vs[i];

		ASSERT_VEC_EQ( 2 * Len, dst, r );

		vec.store_u(dst + 1);

		ASSERT_VEC_EQ( (2 * Len - 1), (dst + 1), r );

	}
};

TN_CASE( sse_svec_store )
{
	svec_store_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_store )
{
	svec_store_test<T, avx_t, N>::run();
}
#endif


template<typename T, typename Kind, index_t Len>
struct svec_add_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs_a[Len];
		for (index_t i = 0; i < Len; ++i) vs_a[i] = T(2 * i + 1);

		T vs_b[Len];
		for (index_t i = 0; i < Len; ++i) vs_b[i] = T(3 * i - 2);

		vec_t va(vs_a);
		vec_t vb(vs_b);

		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		T r[inlen];
		for (index_t i = 0; i < inlen; ++i) r[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r[i] = vs_a[i] + vs_b[i];

		vec_t vr = va + vb;
		const T* pr = (const T*)(&vr);
		ASSERT_VEC_EQ(inlen, pr, r);

		va += vb;
		const T *pa = (const T*)(&va);
		ASSERT_VEC_EQ(inlen, pa, r);
	}
};

TN_CASE( sse_svec_add )
{
	svec_add_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_add )
{
	svec_add_test<T, avx_t, N>::run();
}
#endif


template<typename T, typename Kind, index_t Len>
struct svec_sub_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs_a[Len];
		for (index_t i = 0; i < Len; ++i) vs_a[i] = T(2 * i + 1);

		T vs_b[Len];
		for (index_t i = 0; i < Len; ++i) vs_b[i] = T(3 * i - 2);

		vec_t va(vs_a);
		vec_t vb(vs_b);

		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		T r[inlen];
		for (index_t i = 0; i < inlen; ++i) r[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r[i] = vs_a[i] - vs_b[i];

		vec_t vr = va - vb;
		const T* pr = (const T*)(&vr);
		ASSERT_VEC_EQ(inlen, pr, r);

		va -= vb;
		const T *pa = (const T*)(&va);
		ASSERT_VEC_EQ(inlen, pa, r);
	}
};

TN_CASE( sse_svec_sub )
{
	svec_sub_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_sub )
{
	svec_sub_test<T, avx_t, N>::run();
}
#endif


template<typename T, typename Kind, index_t Len>
struct svec_mul_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs_a[Len];
		for (index_t i = 0; i < Len; ++i) vs_a[i] = T(2 * i + 1);

		T vs_b[Len];
		for (index_t i = 0; i < Len; ++i) vs_b[i] = T(3 * i - 2);

		vec_t va(vs_a);
		vec_t vb(vs_b);

		const index_t inlen = (index_t)(sizeof(vec_t) / sizeof(T));
		T r[inlen];
		for (index_t i = 0; i < inlen; ++i) r[i] = T(0);
		for (index_t i = 0; i < Len; ++i) r[i] = vs_a[i] * vs_b[i];

		vec_t vr = va * vb;
		const T* pr = (const T*)(&vr);
		ASSERT_VEC_EQ(inlen, pr, r);

		va *= vb;
		const T *pa = (const T*)(&va);
		ASSERT_VEC_EQ(inlen, pa, r);
	}
};

TN_CASE( sse_svec_mul )
{
	svec_mul_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_mul )
{
	svec_mul_test<T, avx_t, N>::run();
}
#endif



template<typename T, typename Kind, index_t Len>
struct svec_sum_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs[Len];
		for (index_t i = 0; i < Len; ++i) vs[i] = T(2 * i + 1);
		vec_t vec(vs);

		T s(0);
		for (index_t i = 0; i < Len; ++i) s += vs[i];

		ASSERT_EQ( sum(vec), s );
	}
};

TN_CASE( sse_svec_sum )
{
	svec_sum_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_sum )
{
	svec_sum_test<T, avx_t, N>::run();
}
#endif


template<typename T, typename Kind, index_t Len>
struct svec_dot_test
{
	static void run()
	{
		typedef simd_vec<T, Kind, Len> vec_t;

		T vs_a[Len];
		T vs_b[Len];

		for (index_t i = 0; i < Len; ++i) vs_a[i] = T(2 * i + 1);
		for (index_t i = 0; i < Len; ++i) vs_b[i] = T(3 * i - 2);

		vec_t va(vs_a);
		vec_t vb(vs_b);

		T s(0);
		for (index_t i = 0; i < Len; ++i) s += vs_a[i] * vs_b[i];

		ASSERT_EQ( dot(va, vb), s );
	}
};

TN_CASE( sse_svec_dot )
{
	svec_dot_test<T, sse_t, N>::run();
}

#ifdef LMAT_HAS_AVX
TN_CASE( avx_svec_dot )
{
	svec_dot_test<T, avx_t, N>::run();
}
#endif




AUTO_TPACK( simd_vec_set )
{
	ADD_TN_CASE( sse_svec_set, float,  2 )
	ADD_TN_CASE( sse_svec_set, double, 2 )
	ADD_TN_CASE( sse_svec_set, float,  3 )
	ADD_TN_CASE( sse_svec_set, double, 3 )
	ADD_TN_CASE( sse_svec_set, float,  4 )
	ADD_TN_CASE( sse_svec_set, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_set, float,  2 )
	ADD_TN_CASE( avx_svec_set, double, 2 )
	ADD_TN_CASE( avx_svec_set, float,  3 )
	ADD_TN_CASE( avx_svec_set, double, 3 )
	ADD_TN_CASE( avx_svec_set, float,  4 )
	ADD_TN_CASE( avx_svec_set, double, 4 )
#endif
}

AUTO_TPACK( simd_vec_load )
{
	ADD_TN_CASE( sse_svec_load, float,  2 )
	ADD_TN_CASE( sse_svec_load, double, 2 )
	ADD_TN_CASE( sse_svec_load, float,  3 )
	ADD_TN_CASE( sse_svec_load, double, 3 )
	ADD_TN_CASE( sse_svec_load, float,  4 )
	ADD_TN_CASE( sse_svec_load, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_load, float,  2 )
	ADD_TN_CASE( avx_svec_load, double, 2 )
	ADD_TN_CASE( avx_svec_load, float,  3 )
	ADD_TN_CASE( avx_svec_load, double, 3 )
	ADD_TN_CASE( avx_svec_load, float,  4 )
	ADD_TN_CASE( avx_svec_load, double, 4 )
#endif
}

AUTO_TPACK( simd_vec_store )
{
	ADD_TN_CASE( sse_svec_store, float,  2 )
	ADD_TN_CASE( sse_svec_store, double, 2 )
	ADD_TN_CASE( sse_svec_store, float,  3 )
	ADD_TN_CASE( sse_svec_store, double, 3 )
	ADD_TN_CASE( sse_svec_store, float,  4 )
	ADD_TN_CASE( sse_svec_store, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_store, float,  2 )
	ADD_TN_CASE( avx_svec_store, double, 2 )
	ADD_TN_CASE( avx_svec_store, float,  3 )
	ADD_TN_CASE( avx_svec_store, double, 3 )
	ADD_TN_CASE( avx_svec_store, float,  4 )
	ADD_TN_CASE( avx_svec_store, double, 4 )
#endif
}


AUTO_TPACK( simd_vec_add )
{
	ADD_TN_CASE( sse_svec_add, float,  2 )
	ADD_TN_CASE( sse_svec_add, double, 2 )
	ADD_TN_CASE( sse_svec_add, float,  3 )
	ADD_TN_CASE( sse_svec_add, double, 3 )
	ADD_TN_CASE( sse_svec_add, float,  4 )
	ADD_TN_CASE( sse_svec_add, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_add, float,  2 )
	ADD_TN_CASE( avx_svec_add, double, 2 )
	ADD_TN_CASE( avx_svec_add, float,  3 )
	ADD_TN_CASE( avx_svec_add, double, 3 )
	ADD_TN_CASE( avx_svec_add, float,  4 )
	ADD_TN_CASE( avx_svec_add, double, 4 )
#endif
}


AUTO_TPACK( simd_vec_sub )
{
	ADD_TN_CASE( sse_svec_sub, float,  2 )
	ADD_TN_CASE( sse_svec_sub, double, 2 )
	ADD_TN_CASE( sse_svec_sub, float,  3 )
	ADD_TN_CASE( sse_svec_sub, double, 3 )
	ADD_TN_CASE( sse_svec_sub, float,  4 )
	ADD_TN_CASE( sse_svec_sub, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_sub, float,  2 )
	ADD_TN_CASE( avx_svec_sub, double, 2 )
	ADD_TN_CASE( avx_svec_sub, float,  3 )
	ADD_TN_CASE( avx_svec_sub, double, 3 )
	ADD_TN_CASE( avx_svec_sub, float,  4 )
	ADD_TN_CASE( avx_svec_sub, double, 4 )
#endif
}


AUTO_TPACK( simd_vec_mul )
{
	ADD_TN_CASE( sse_svec_mul, float,  2 )
	ADD_TN_CASE( sse_svec_mul, double, 2 )
	ADD_TN_CASE( sse_svec_mul, float,  3 )
	ADD_TN_CASE( sse_svec_mul, double, 3 )
	ADD_TN_CASE( sse_svec_mul, float,  4 )
	ADD_TN_CASE( sse_svec_mul, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_mul, float,  2 )
	ADD_TN_CASE( avx_svec_mul, double, 2 )
	ADD_TN_CASE( avx_svec_mul, float,  3 )
	ADD_TN_CASE( avx_svec_mul, double, 3 )
	ADD_TN_CASE( avx_svec_mul, float,  4 )
	ADD_TN_CASE( avx_svec_mul, double, 4 )
#endif
}


AUTO_TPACK( simd_vec_sum )
{
	ADD_TN_CASE( sse_svec_sum, float,  2 )
	ADD_TN_CASE( sse_svec_sum, double, 2 )
	ADD_TN_CASE( sse_svec_sum, float,  3 )
	ADD_TN_CASE( sse_svec_sum, double, 3 )
	ADD_TN_CASE( sse_svec_sum, float,  4 )
	ADD_TN_CASE( sse_svec_sum, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_sum, float,  2 )
	ADD_TN_CASE( avx_svec_sum, double, 2 )
	ADD_TN_CASE( avx_svec_sum, float,  3 )
	ADD_TN_CASE( avx_svec_sum, double, 3 )
	ADD_TN_CASE( avx_svec_sum, float,  4 )
	ADD_TN_CASE( avx_svec_sum, double, 4 )
#endif
}


AUTO_TPACK( simd_vec_dot )
{
	ADD_TN_CASE( sse_svec_dot, float,  2 )
	ADD_TN_CASE( sse_svec_dot, double, 2 )
	ADD_TN_CASE( sse_svec_dot, float,  3 )
	ADD_TN_CASE( sse_svec_dot, double, 3 )
	ADD_TN_CASE( sse_svec_dot, float,  4 )
	ADD_TN_CASE( sse_svec_dot, double, 4 )

#ifdef LMAT_HAS_AVX
	ADD_TN_CASE( avx_svec_dot, float,  2 )
	ADD_TN_CASE( avx_svec_dot, double, 2 )
	ADD_TN_CASE( avx_svec_dot, float,  3 )
	ADD_TN_CASE( avx_svec_dot, double, 3 )
	ADD_TN_CASE( avx_svec_dot, float,  4 )
	ADD_TN_CASE( avx_svec_dot, double, 4 )
#endif
}

