/**
 * @file const_matrix.h
 *
 * Constant matrix class
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_CONST_MATRIX_H_
#define LIGHTMAT_CONST_MATRIX_H_

#include <light_mat/matrix/matrix_fill.h>
#include <light_mat/matrix/matrix_shape.h>

namespace lmat
{

	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<const_matrix<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = true;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols>
	struct is_linear_accessible<const_matrix<T, CTRows, CTCols> >
	{
		static const bool value = true;
	};

	struct matrix_fill_policy { };

	template<typename T, int CTRows, int CTCols, class DMat>
	struct default_matrix_eval_policy<const_matrix<T, CTRows, CTCols>, DMat>
	{
		typedef matrix_fill_policy type;
	};


	template<typename T, int CTRows, int CTCols>
	class const_matrix : public IMatrixView<const_matrix<T, CTRows, CTCols>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

		LMAT_ENSURE_INLINE const_matrix(const index_type m, const index_type n, const T& val)
		: m_shape(m, n), m_val(val) { }

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE T value() const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE T elem(const index_type, const index_type) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE T operator[] (const index_type) const
		{
			return m_val;
		}

	private:
		matrix_shape<CTRows, CTCols> m_shape;
		const T m_val;
	};


	// Evaluation

	template<typename T, int M, int N, class Dst>
	LMAT_ENSURE_INLINE
	void evaluate(const const_matrix<T, M, N>& mat, IDenseMatrix<Dst, T>& dst, matrix_fill_policy)
	{
		fill(dst, mat.value());
	}


	// Transpose

	template<typename T, int M, int N, typename Arg_HP>
	struct unary_expr_map<transpose_t, Arg_HP, const_matrix<T, M, N> >
	{
		typedef const_matrix<T, N, M> type;

		LMAT_ENSURE_INLINE
		static type get(transpose_t, const arg_forwarder<Arg_HP, const_matrix<T, M, N> >& arg_fwd)
		{
			const const_matrix<T, M, N>& arg = arg_fwd.arg;
			return type(arg.ncolumns(), arg.nrows(), arg.value());
		}
	};

}


#endif 
