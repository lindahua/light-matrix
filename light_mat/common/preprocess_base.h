/**
 * @file preprocess_base.h
 *
 * @brief Useful preprocessors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_PREPROCESS_BASE_H_
#define LIGHTMAT_PREPROCESS_BASE_H_

#include <light_mat/common/prim_types.h>

#define LMAT_REPEAT_ARGS_1(A) A
#define LMAT_REPEAT_ARGS_2(A) A, A
#define LMAT_REPEAT_ARGS_3(A) A, A, A
#define LMAT_REPEAT_ARGS_4(A) A, A, A, A
#define LMAT_REPEAT_ARGS_5(A) A, A, A, A, A
#define LMAT_REPEAT_ARGS_6(A) A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_7(A) A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_8(A) A, A, A, A, A, A, A, A

#define LMAT_REPEAT_ARGS_P1(A) A##1
#define LMAT_REPEAT_ARGS_P2(A) A##1, A##2
#define LMAT_REPEAT_ARGS_P3(A) A##1, A##2, A##3
#define LMAT_REPEAT_ARGS_P4(A) A##1, A##2, A##3, A##4
#define LMAT_REPEAT_ARGS_P5(A) A##1, A##2, A##3, A##4, A##5
#define LMAT_REPEAT_ARGS_P6(A) A##1, A##2, A##3, A##4, A##5, A##6
#define LMAT_REPEAT_ARGS_P7(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7
#define LMAT_REPEAT_ARGS_P8(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8

#define LMAT_REPEAT_ARGS_F1(F)
#define LMAT_REPEAT_ARGS_F2(F) F(1), F(2)
#define LMAT_REPEAT_ARGS_F3(F) F(1), F(2), F(3)
#define LMAT_REPEAT_ARGS_F4(F) F(1), F(2), F(3), F(4)
#define LMAT_REPEAT_ARGS_F5(F) F(1), F(2), F(3), F(4), F(5)
#define LMAT_REPEAT_ARGS_F6(F) F(1), F(2), F(3), F(4), F(5), F(6)
#define LMAT_REPEAT_ARGS_F7(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7)
#define LMAT_REPEAT_ARGS_F8(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8)

#define LMAT_REPEAT_STATEMENTS_F1(S) S(1);
#define LMAT_REPEAT_STATEMENTS_F2(S) \
		S(1); \
		S(2);
#define LMAT_REPEAT_STATEMENTS_F3(S) \
		S(1); \
		S(2); \
		S(3);
#define LMAT_REPEAT_STATEMENTS_F4(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4);
#define LMAT_REPEAT_STATEMENTS_F5(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5);
#define LMAT_REPEAT_STATEMENTS_F6(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6);
#define LMAT_REPEAT_STATEMENTS_F7(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7);
#define LMAT_REPEAT_STATEMENTS_F8(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8);

#endif
