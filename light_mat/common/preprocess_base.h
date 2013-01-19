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
#define LMAT_REPEAT_ARGS_9(A) A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_10(A) A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_11(A) A, A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_12(A) A, A, A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_13(A) A, A, A, A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_14(A) A, A, A, A, A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_15(A) A, A, A, A, A, A, A, A, A, A, A, A, A, A, A
#define LMAT_REPEAT_ARGS_16(A) A, A, A, A, A, A, A, A, A, A, A, A, A, A, A, A

#define LMAT_REPEAT_ARGS_P1(A) A##1
#define LMAT_REPEAT_ARGS_P2(A) A##1, A##2
#define LMAT_REPEAT_ARGS_P3(A) A##1, A##2, A##3
#define LMAT_REPEAT_ARGS_P4(A) A##1, A##2, A##3, A##4
#define LMAT_REPEAT_ARGS_P5(A) A##1, A##2, A##3, A##4, A##5
#define LMAT_REPEAT_ARGS_P6(A) A##1, A##2, A##3, A##4, A##5, A##6
#define LMAT_REPEAT_ARGS_P7(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7
#define LMAT_REPEAT_ARGS_P8(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8
#define LMAT_REPEAT_ARGS_P9(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9
#define LMAT_REPEAT_ARGS_P10(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10
#define LMAT_REPEAT_ARGS_P11(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11
#define LMAT_REPEAT_ARGS_P12(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11, A##12
#define LMAT_REPEAT_ARGS_P13(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11, A##12, A##13
#define LMAT_REPEAT_ARGS_P14(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11, A##12, A##13, A##14
#define LMAT_REPEAT_ARGS_P15(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11, A##12, A##13, A##14, A##15
#define LMAT_REPEAT_ARGS_P16(A) A##1, A##2, A##3, A##4, A##5, A##6, A##7, A##8, A##9, A##10, A##11, A##12, A##13, A##14, A##15, A##16

#define LMAT_REPEAT_ARGS_F1(F)
#define LMAT_REPEAT_ARGS_F2(F) F(1), F(2)
#define LMAT_REPEAT_ARGS_F3(F) F(1), F(2), F(3)
#define LMAT_REPEAT_ARGS_F4(F) F(1), F(2), F(3), F(4)
#define LMAT_REPEAT_ARGS_F5(F) F(1), F(2), F(3), F(4), F(5)
#define LMAT_REPEAT_ARGS_F6(F) F(1), F(2), F(3), F(4), F(5), F(6)
#define LMAT_REPEAT_ARGS_F7(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7)
#define LMAT_REPEAT_ARGS_F8(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8)
#define LMAT_REPEAT_ARGS_F9(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9)
#define LMAT_REPEAT_ARGS_F10(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10)
#define LMAT_REPEAT_ARGS_F11(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11)
#define LMAT_REPEAT_ARGS_F12(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11), F(12)
#define LMAT_REPEAT_ARGS_F13(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11), F(12), F(13)
#define LMAT_REPEAT_ARGS_F14(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11), F(12), F(13), F(14)
#define LMAT_REPEAT_ARGS_F15(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11), F(12), F(13), F(14), F(15)
#define LMAT_REPEAT_ARGS_F16(F) F(1), F(2), F(3), F(4), F(5), F(6), F(7), F(8), F(9), F(10), F(11), F(12), F(13), F(14), F(15), F(16)

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
#define LMAT_REPEAT_STATEMENTS_F9(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9);
#define LMAT_REPEAT_STATEMENTS_F10(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10);
#define LMAT_REPEAT_STATEMENTS_F11(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11);
#define LMAT_REPEAT_STATEMENTS_F12(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11); \
		S(12);
#define LMAT_REPEAT_STATEMENTS_F13(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11); \
		S(12); \
		S(13);
#define LMAT_REPEAT_STATEMENTS_F14(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11); \
		S(12); \
		S(13); \
		S(14);
#define LMAT_REPEAT_STATEMENTS_F15(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11); \
		S(12); \
		S(13); \
		S(14); \
		S(15);
#define LMAT_REPEAT_STATEMENTS_F16(S) \
		S(1); \
		S(2); \
		S(3); \
		S(4); \
		S(5); \
		S(6); \
		S(7); \
		S(8); \
		S(9); \
		S(10); \
		S(11); \
		S(12); \
		S(13); \
		S(14); \
		S(15); \
		S(16);

#endif
