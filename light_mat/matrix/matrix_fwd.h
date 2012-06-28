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

	// forward declaration of some important types

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim> class dense_matrix; // guranteed to be aligned
	template<typename T, int CTRows=DynamicDim> class dense_col;
	template<typename T, int CTCols=DynamicDim> class dense_row;

	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim> class cref_matrix;
	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim> class ref_matrix;

	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim> class aligned_cref_matrix;
	template<typename T, int RowDim=DynamicDim, int ColDim=DynamicDim> class aligned_ref_matrix;

	template<typename T, int CTRows=DynamicDim> class cref_col;
	template<typename T, int CTRows=DynamicDim> class ref_col;
	template<typename T, int CTCols=DynamicDim> class cref_row;
	template<typename T, int CTCols=DynamicDim> class ref_row;

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim> class cref_matrix_ex;
	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim> class ref_matrix_ex;

	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim> class cref_grid2d;
	template<typename T, int CTRows=DynamicDim, int CTCols=DynamicDim> class ref_grid2d;

	// auxiliary structures

	template<class Mat, typename RowRange> struct colviews;
	template<class Mat, typename ColRange> struct rowviews;
	template<class Mat, typename RowRange, typename ColRange> struct subviews;

	template<class Mat> struct transposed;
	template<class Mat> class transpose_expr;

	template<class Mat> struct mat_evaluator;

	// forward declaration of useful memory operations

	template<typename T, class LMat, class RMat>
	inline bool is_equal(const IMatrixView<LMat, T>& A, const IMatrixView<RMat, T>& B);

	template<typename T, class LMat, class RMat>
	inline bool is_approx(const IMatrixView<LMat, T>& A, const IMatrixView<RMat, T>& B, const T& tol);

	template<typename T, class SExpr, class DMat>
	inline void evaluate_to(const IMatrixXpr<SExpr, T>& src, IDenseMatrix<DMat, T>& dst);

	template<typename T, class SMat, class DMat>
	inline void copy(const IMatrixView<SMat, T>& src, IDenseMatrix<DMat, T>& dst);

	template<typename T, class DMat>
	inline void fill(IDenseMatrix<DMat, T>& dst, const T& v);

	template<typename T, class DMat>
	inline void zero(IDenseMatrix<DMat, T>& dst);
}


#endif
