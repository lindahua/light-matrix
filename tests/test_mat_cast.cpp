/**
 * @file test_dense_cast.cpp
 *
 * Unit testing of casting for dense matrices
 *
 * @author Dahua Lin
 */

#include "test_base.h"

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/matrix_cast.h>

using namespace lmat;
using namespace lmat::test;

const index_t DM = 5;
const index_t DN = 6;


#define TEST_MAT_CAST(S, T, to_fun) \
	MN_CASE( mat_cast, S##_to_##T ) { \
		const index_t m = M == 0 ? DM : M; \
		const index_t n = N == 0 ? DN : N; \
		dense_matrix<S, M, N> smat(m, n); \
		for (index_t i = 0; i < m * n; ++i) { \
			smat[i] = S(i + 1); \
		} \
		dense_matrix<T, M, N> rmat(m, n); \
		for (index_t i = 0; i < m * n; ++i) { \
			rmat[i] = static_cast<T>(smat[i]); \
		} \
		dense_matrix<T, M, N> tmat1 = cast(smat, type<T>()); \
		ASSERT_EQ( tmat1.nrows(), m ); \
		ASSERT_EQ( tmat1.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, tmat1, rmat ); \
		dense_matrix<T, M, N> tmat2 = to_fun(smat); \
		ASSERT_EQ( tmat2.nrows(), m ); \
		ASSERT_EQ( tmat2.ncolumns(), n ); \
		ASSERT_MAT_EQ( m, n, tmat2, rmat ); \
	} \
	BEGIN_TPACK( mat_cast_##S##_to_##T ) \
		ADD_MN_CASE_3X3( mat_cast, S##_to_##T, DM, DN ) \
	END_TPACK


#define ADD_CAST_TPACK( S, T ) ADD_TPACK( mat_cast_##S##_to_##T )


TEST_MAT_CAST( int8_t, double, to_f64 )
TEST_MAT_CAST( uint8_t, double, to_f64 )
TEST_MAT_CAST( int16_t, double, to_f64 )
TEST_MAT_CAST( uint16_t, double, to_f64 )
TEST_MAT_CAST( int32_t, double, to_f64 )
TEST_MAT_CAST( uint32_t, double, to_f64 )
TEST_MAT_CAST( float, double, to_f64 )
TEST_MAT_CAST( double, double, to_f64 )

TEST_MAT_CAST( int8_t, float, to_f32 )
TEST_MAT_CAST( uint8_t, float, to_f32 )
TEST_MAT_CAST( int16_t, float, to_f32 )
TEST_MAT_CAST( uint16_t, float, to_f32 )
TEST_MAT_CAST( int32_t, float, to_f32 )
TEST_MAT_CAST( uint32_t, float, to_f32 )
TEST_MAT_CAST( float, float, to_f32 )
TEST_MAT_CAST( double, float, to_f32 )

TEST_MAT_CAST( int8_t, int32_t, to_i32 )
TEST_MAT_CAST( uint8_t, int32_t, to_i32 )
TEST_MAT_CAST( int16_t, int32_t, to_i32 )
TEST_MAT_CAST( uint16_t, int32_t, to_i32 )
TEST_MAT_CAST( int32_t, int32_t, to_i32 )
TEST_MAT_CAST( uint32_t, int32_t, to_i32 )
TEST_MAT_CAST( float, int32_t, to_i32 )
TEST_MAT_CAST( double, int32_t, to_i32 )

TEST_MAT_CAST( int8_t, uint32_t, to_u32 )
TEST_MAT_CAST( uint8_t, uint32_t, to_u32 )
TEST_MAT_CAST( int16_t, uint32_t, to_u32 )
TEST_MAT_CAST( uint16_t, uint32_t, to_u32 )
TEST_MAT_CAST( int32_t, uint32_t, to_u32 )
TEST_MAT_CAST( uint32_t, uint32_t, to_u32 )
TEST_MAT_CAST( float, uint32_t, to_u32 )
TEST_MAT_CAST( double, uint32_t, to_u32 )

TEST_MAT_CAST( int8_t, int16_t, to_i16 )
TEST_MAT_CAST( uint8_t, int16_t, to_i16 )
TEST_MAT_CAST( int16_t, int16_t, to_i16 )
TEST_MAT_CAST( uint16_t, int16_t, to_i16 )
TEST_MAT_CAST( int32_t, int16_t, to_i16 )
TEST_MAT_CAST( uint32_t, int16_t, to_i16 )
TEST_MAT_CAST( float, int16_t, to_i16 )
TEST_MAT_CAST( double, int16_t, to_i16 )

TEST_MAT_CAST( int8_t, uint16_t, to_u16 )
TEST_MAT_CAST( uint8_t, uint16_t, to_u16 )
TEST_MAT_CAST( int16_t, uint16_t, to_u16 )
TEST_MAT_CAST( uint16_t, uint16_t, to_u16 )
TEST_MAT_CAST( int32_t, uint16_t, to_u16 )
TEST_MAT_CAST( uint32_t, uint16_t, to_u16 )
TEST_MAT_CAST( float, uint16_t, to_u16 )
TEST_MAT_CAST( double, uint16_t, to_u16 )

TEST_MAT_CAST( int8_t, int8_t, to_i8 )
TEST_MAT_CAST( uint8_t, int8_t, to_i8 )
TEST_MAT_CAST( int16_t, int8_t, to_i8 )
TEST_MAT_CAST( uint16_t, int8_t, to_i8 )
TEST_MAT_CAST( int32_t, int8_t, to_i8 )
TEST_MAT_CAST( uint32_t, int8_t, to_i8 )
TEST_MAT_CAST( float, int8_t, to_i8 )
TEST_MAT_CAST( double, int8_t, to_i8 )

TEST_MAT_CAST( int8_t, uint8_t, to_u8 )
TEST_MAT_CAST( uint8_t, uint8_t, to_u8 )
TEST_MAT_CAST( int16_t, uint8_t, to_u8 )
TEST_MAT_CAST( uint16_t, uint8_t, to_u8 )
TEST_MAT_CAST( int32_t, uint8_t, to_u8 )
TEST_MAT_CAST( uint32_t, uint8_t, to_u8 )
TEST_MAT_CAST( float, uint8_t, to_u8 )
TEST_MAT_CAST( double, uint8_t, to_u8 )


BEGIN_MAIN_SUITE
	ADD_CAST_TPACK( int8_t, double )
	ADD_CAST_TPACK( uint8_t, double )
	ADD_CAST_TPACK( int16_t, double )
	ADD_CAST_TPACK( uint16_t, double )
	ADD_CAST_TPACK( int32_t, double )
	ADD_CAST_TPACK( uint32_t, double )
	ADD_CAST_TPACK( float, double )
	ADD_CAST_TPACK( double, double )

	ADD_CAST_TPACK( int8_t, float )
	ADD_CAST_TPACK( uint8_t, float )
	ADD_CAST_TPACK( int16_t, float )
	ADD_CAST_TPACK( uint16_t, float )
	ADD_CAST_TPACK( int32_t, float )
	ADD_CAST_TPACK( uint32_t, float )
	ADD_CAST_TPACK( float, float )
	ADD_CAST_TPACK( double, float )

	ADD_CAST_TPACK( int8_t, int32_t )
	ADD_CAST_TPACK( uint8_t, int32_t )
	ADD_CAST_TPACK( int16_t, int32_t )
	ADD_CAST_TPACK( uint16_t, int32_t )
	ADD_CAST_TPACK( int32_t, int32_t )
	ADD_CAST_TPACK( uint32_t, int32_t )
	ADD_CAST_TPACK( float, int32_t )
	ADD_CAST_TPACK( double, int32_t)

	ADD_CAST_TPACK( int8_t, uint32_t )
	ADD_CAST_TPACK( uint8_t, uint32_t )
	ADD_CAST_TPACK( int16_t, uint32_t )
	ADD_CAST_TPACK( uint16_t, uint32_t )
	ADD_CAST_TPACK( int32_t, uint32_t )
	ADD_CAST_TPACK( uint32_t, uint32_t )
	ADD_CAST_TPACK( float, uint32_t )
	ADD_CAST_TPACK( double, uint32_t )

	ADD_CAST_TPACK( int8_t, int16_t )
	ADD_CAST_TPACK( uint8_t, int16_t )
	ADD_CAST_TPACK( int16_t, int16_t )
	ADD_CAST_TPACK( uint16_t, int16_t )
	ADD_CAST_TPACK( int32_t, int16_t )
	ADD_CAST_TPACK( uint32_t, int16_t )
	ADD_CAST_TPACK( float, int16_t )
	ADD_CAST_TPACK( double, int16_t )

	ADD_CAST_TPACK( int8_t, uint16_t )
	ADD_CAST_TPACK( uint8_t, uint16_t )
	ADD_CAST_TPACK( int16_t, uint16_t )
	ADD_CAST_TPACK( uint16_t, uint16_t )
	ADD_CAST_TPACK( int32_t, uint16_t )
	ADD_CAST_TPACK( uint32_t, uint16_t )
	ADD_CAST_TPACK( float, uint16_t )
	ADD_CAST_TPACK( double, uint16_t )

	ADD_CAST_TPACK( int8_t, int8_t )
	ADD_CAST_TPACK( uint8_t, int8_t )
	ADD_CAST_TPACK( int16_t, int8_t )
	ADD_CAST_TPACK( uint16_t, int8_t )
	ADD_CAST_TPACK( int32_t, int8_t )
	ADD_CAST_TPACK( uint32_t, int8_t )
	ADD_CAST_TPACK( float, int8_t )
	ADD_CAST_TPACK( double, int8_t )

	ADD_CAST_TPACK( int8_t, uint8_t )
	ADD_CAST_TPACK( uint8_t, uint8_t )
	ADD_CAST_TPACK( int16_t, uint8_t )
	ADD_CAST_TPACK( uint16_t, uint8_t )
	ADD_CAST_TPACK( int32_t, uint8_t )
	ADD_CAST_TPACK( uint32_t, uint8_t )
	ADD_CAST_TPACK( float, uint8_t )
	ADD_CAST_TPACK( double, uint8_t )
END_MAIN_SUITE





