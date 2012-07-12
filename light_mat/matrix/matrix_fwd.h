/**
 * @file matrix_fwd.h
 *
 * Forward declaration of Matrix-related classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_FWD_H_
#define LIGHTMAT_MATRIX_FWD_H_

#include <light_mat/core/basic_defs.h>
#include <light_mat/core/range.h>

namespace lmat
{
	const int DynamicDim = 0;

	template<int N>
	struct fixed_dim { static const int value = N; };

	// forward declaration of concepts


	/****************************************************************
	 *
	 *   A specialized version of matrix_traits<C> must be provided,
	 *   which should contain the following static members:
	 *
	 *   - num_dimensions:	an int value, which must be set to 2
	 *   					(reserved for future extension)
	 *
	 *   - compile_time_num_rows:	compile-time number of rows
	 *   - compile_time_num_cols:	compile-time number of columns
	 *   - is_readonly:				whether the contents can be modified
	 *
	 *	 - value_type:			the type of element value
	 *
	 ****************************************************************/

	template<class Derived> struct matrix_traits;

	template<class Derived, typename T> class IMatrixXpr;
	template<class Derived, typename T> class IMatrixView;
	template<class Derived, typename T> class IDenseMatrix;

	// alignment tags

	struct unaligned { };
	struct base_aligned { };
	struct percol_aligned { };

	template<typename T>
	struct is_align_tag { static const bool value = false; };

	template<> struct is_align_tag<unaligned> { static const bool value = true; };
	template<> struct is_align_tag<base_aligned> { static const bool value = true; };
	template<> struct is_align_tag<percol_aligned> { static const bool value = true; };

	// forward declaration of some important types

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim, typename Align=base_aligned>
	class dense_matrix;

	template<typename T, int CTRows=DynamicDim, typename Align=base_aligned> class dense_col;
	template<typename T, int CTCols=DynamicDim, typename Align=base_aligned> class dense_row;

	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim, typename Align=unaligned>
	class cref_matrix;

	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim, typename Align=unaligned>
	class ref_matrix;

	template<typename T, int CTRows=DynamicDim, typename Align=unaligned> class cref_col;
	template<typename T, int CTRows=DynamicDim, typename Align=unaligned> class ref_col;
	template<typename T, int CTCols=DynamicDim, typename Align=unaligned> class cref_row;
	template<typename T, int CTCols=DynamicDim, typename Align=unaligned> class ref_row;

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim, typename Align=unaligned>
	class cref_matrix_ex;

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim, typename Align=unaligned>
	class ref_matrix_ex;

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim>
	class const_matrix;

	// auxiliary structures

	template<class Mat, typename RowRange> struct colviews;
	template<class Mat, typename ColRange> struct rowviews;
	template<class Mat, typename RowRange, typename ColRange> struct subviews;

	// expressions

	template<class Fun, class Arg> class unary_ewise_expr;
	template<class Fun, class Arg1, class Arg2> class binary_ewise_expr;
	template<class Fun, class Arg2> class binary_fix1_ewise_expr;
	template<class Fun, class Arg1> class binary_fix2_ewise_expr;

	template<class Col, int N=DynamicDim> class repeat_col_expr;
	template<class Row, int M=DynamicDim> class repeat_row_expr;

	struct rowwise { };
	struct colwise { };

	template<class Fun, class Arg> class colwise_reduce_expr;
	template<class Fun, class Arg> class rowwise_reduce_expr;

	template<class Mat> struct transposed;
	template<class Mat> class transpose_expr;

	// evaluation

	template<class Expr> struct linear_vector_evaluator;
	template<class Expr> struct percol_vector_evaluator;
	template<class Expr> struct linear_simd_vector_evaluator;
	template<class Expr> struct percol_simd_vector_evaluator;

	template<typename T, int CTSize> class linear_eval_context;
	template<typename T, int CTRows, int CTCols> class percol_eval_context;

	template<class Expr> struct linear_eval_cost;
	template<class Expr> struct percol_eval_cost;

	template<typename T, class SExpr, class DMat>
	inline void evaluate_to(const IMatrixXpr<SExpr, T>& src, IDenseMatrix<DMat, T>& dst);

	template<class SMat, typename T> struct cast_expr_map;

}

// Useful macros

#define LMAT_MATRIX_TYPEDEFS0(TName, prefix) \
	typedef TName<double>   prefix##_f64; \
	typedef TName<float>    prefix##_f32; \
	typedef TName<int32_t>  prefix##_i32; \
	typedef TName<uint32_t> prefix##_u32; \
	typedef TName<int16_t>  prefix##_i16; \
	typedef TName<uint16_t> prefix##_u16; \
	typedef TName<int8_t>   prefix##_i8; \
	typedef TName<uint8_t>  prefix##_u8; \
	typedef TName<bool>     prefix##_bool;

#define LMAT_MATRIX_TYPEDEFS1(TName, prefix, Dim) \
	typedef TName<double,   Dim> prefix##_f64; \
	typedef TName<float,    Dim> prefix##_f32; \
	typedef TName<int32_t,  Dim> prefix##_i32; \
	typedef TName<uint32_t, Dim> prefix##_u32; \
	typedef TName<int16_t,  Dim> prefix##_i16; \
	typedef TName<uint16_t, Dim> prefix##_u16; \
	typedef TName<int8_t,   Dim> prefix##_i8; \
	typedef TName<uint8_t,  Dim> prefix##_u8; \
	typedef TName<bool,     Dim> prefix##_bool;

#define LMAT_MATRIX_TYPEDEFS2(TName, prefix, RDim, CDim) \
	typedef TName<double,   RDim, CDim> prefix##_f64; \
	typedef TName<float,    RDim, CDim> prefix##_f32; \
	typedef TName<int32_t,  RDim, CDim> prefix##_i32; \
	typedef TName<uint32_t, RDim, CDim> prefix##_u32; \
	typedef TName<int16_t,  RDim, CDim> prefix##_i16; \
	typedef TName<uint16_t, RDim, CDim> prefix##_u16; \
	typedef TName<int8_t,   RDim, CDim> prefix##_i8; \
	typedef TName<uint8_t,  RDim, CDim> prefix##_u8; \
	typedef TName<bool,     RDim, CDim> prefix##_bool;


#endif


