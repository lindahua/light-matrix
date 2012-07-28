/**
 * @file repeat_vecs.h
 *
 * Matrices formed by repeating rows/columns
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REPEAT_VECS_H_
#define LIGHTMAT_REPEAT_VECS_H_

#include <light_mat/matrix/matrix_expr_base.h>
#include <light_mat/matrix/matrix_shape.h>

namespace lmat
{
	/********************************************
	 *
	 *  Expression type mapping
	 *
	 ********************************************/

	template<class Col, int N>
	struct repcol_type_map
	{
		typedef repeat_col_expr<Col, N> type;
	};

	template<class Row, int M>
	struct reprow_type_map
	{
		typedef repeat_row_expr<Row, M> type;
	};

	template<typename T, int Mc, int Nc, int N>
	struct repcol_type_map<const_matrix<T, Mc, Nc>, N>
	{
		typedef const_matrix<T, Mc, N> type;
	};

	template<typename T, int Mc, int Nc, int M>
	struct reprow_type_map<const_matrix<T, Mc, Nc>, M>
	{
		typedef const_matrix<T, M, Nc> type;
	};


	/********************************************
	 *
	 *  Expression traits
	 *
	 ********************************************/

	template<class Col, int N>
	struct matrix_traits<repeat_col_expr<Col, N> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Col>::value;
		static const int compile_time_num_cols = N;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Col>::value_type value_type;
	};

	template<class Row, int M>
	struct matrix_traits<repeat_row_expr<Row, M> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = M;
		static const int compile_time_num_cols = ct_cols<Row>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Row>::value_type value_type;
	};

	template<class Col, int N, class Dst>
	struct default_evalctx<repeat_col_expr<Col, N>, Dst>
	{
		static const int N_ = binary_ctdim<N, ct_cols<Dst>::value>::value;
		typedef repcols_evalctx<Col, N_, Dst> type;
	};

	template<class Row, int M, class Dst>
	struct default_evalctx<repeat_row_expr<Row, M>, Dst>
	{
		static const int M_ = binary_ctdim<M, ct_rows<Dst>::value>::value;
		typedef reprows_evalctx<Row, M_, Dst> type;
	};


	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	template<class Col, int N>
	class repeat_col_expr
	: public IMatrixXpr<repeat_col_expr<Col, N>, typename matrix_traits<Col>::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Col>::value, "Col must be a matrix expression class.");
#endif

		typedef typename unwrapped_expr<Col>::type arg_expr_t;

	public:
		typedef typename matrix_traits<Col>::value_type value_type;

		LMAT_ENSURE_INLINE
		repeat_col_expr(const Col& col, const index_t n)
		: m_col(col), m_ncols(n)
		{
			check_arg( is_column(col), "repeat_col: the input col is NOT a column." );
		}

		LMAT_ENSURE_INLINE const arg_expr_t& column() const
		{
			return m_col.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * ncolumns();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return column().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols.dim();
		}

	private:
		obj_wrapper<Col> m_col;
		single_dim<N> m_ncols;
	};


	// repeat_row_expr

	template<class Row, int M>
	class repeat_row_expr
	: public IMatrixXpr<repeat_row_expr<Row, M>, typename matrix_traits<Row>::value_type>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Row>::value, "Row must be a matrix expression class.");
#endif

		typedef typename unwrapped_expr<Row>::type arg_expr_t;

	public:
		typedef typename matrix_traits<Row>::value_type value_type;

		LMAT_ENSURE_INLINE
		repeat_row_expr(const Row& row, const index_t m)
		: m_row(row), m_nrows(m)
		{
			check_arg( is_row(row), "repeat_row: the input row is NOT a row." );
		}

		LMAT_ENSURE_INLINE const arg_expr_t& row() const
		{
			return m_row.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * ncolumns();
		}

		LMAT_ENSURE_INLINE size_t size() const
		{
			return static_cast<size_t>(nelems());
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_nrows.dim();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return row().ncolumns();
		}

	private:
		obj_wrapper<Row> m_row;
		single_dim<M> m_nrows;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename T, class Col>
	LMAT_ENSURE_INLINE
	inline typename repcol_type_map<Col, DynamicDim>::type
	repcol(const IMatrixXpr<Col, T>& col, const index_t n)
	{
		return repeat_col_expr<Col, DynamicDim>(col.derived(), n);
	}

	template<typename T, class Col, int N>
	LMAT_ENSURE_INLINE
	inline typename repcol_type_map<Col, N>::type
	repcol(const IMatrixXpr<Col, T>& col, fixed_dim<N>)
	{
		return repeat_col_expr<Col, N>(col.derived(), N);
	}

	template<typename T, class Row>
	LMAT_ENSURE_INLINE
	inline typename reprow_type_map<Row, DynamicDim>::type
	reprow(const IMatrixXpr<Row, T>& row, const index_t m)
	{
		return repeat_row_expr<Row, DynamicDim>(row.derived(), m);
	}

	template<typename T, class Row, int M>
	LMAT_ENSURE_INLINE
	inline typename reprow_type_map<Row, M>::type
	reprow(const IMatrixXpr<Row, T>& row, fixed_dim<M>)
	{
		return repeat_row_expr<Row, M>(row.derived(), M);
	}


	template<typename T, int Mc, int Nc>
	LMAT_ENSURE_INLINE
	inline const_matrix<T, Mc, DynamicDim>
	repcol(const const_matrix<T, Mc, Nc>& col, const index_t n)
	{
		return const_matrix<T, Mc, DynamicDim>(col.nrows(), n, col.value());
	}

	template<typename T, int Mc, int Nc, int N>
	LMAT_ENSURE_INLINE
	inline const_matrix<T, Mc, N>
	repcol(const const_matrix<T, Mc, Nc>& col, fixed_dim<N>)
	{
		return const_matrix<T, Mc, N>(col.nrows(), N, col.value());
	}

	template<typename T, int Mc, int Nc>
	LMAT_ENSURE_INLINE
	inline const_matrix<T, DynamicDim, Nc>
	reprow(const const_matrix<T, Mc, Nc>& row, const index_t m)
	{
		return const_matrix<T, DynamicDim, Nc>(m, row.ncolumns(), row.value());
	}

	template<typename T, int Mc, int Nc, int M>
	LMAT_ENSURE_INLINE
	inline const_matrix<T, M, Nc>
	reprow(const const_matrix<T, Mc, Nc>& row, fixed_dim<M>)
	{
		return const_matrix<T, M, Nc>(M, row.ncolumns(), row.value());
	}


	/********************************************
	 *
	 *  Transpose
	 *
	 ********************************************/

	template<class Col, int N>
	struct transpose_expr_map<repeat_col_expr<Col, N> >
	{
		typedef repeat_row_expr<
				embed_mat<typename transpose_expr_map<Col>::type>, N> type;

		LMAT_ENSURE_INLINE
		static type get(const repeat_col_expr<Col, N>& expr)
		{
			return type(embed(expr.column().trans()), expr.ncolumns());
		}
	};


	template<class Row, int M>
	struct transpose_expr_map<repeat_row_expr<Row, M> >
	{
		typedef repeat_col_expr<
				embed_mat<typename transpose_expr_map<Row>::type>, M> type;

		LMAT_ENSURE_INLINE
		static type get(const repeat_row_expr<Row, M>& expr)
		{
			return type(embed(expr.row().trans()), expr.nrows());
		}
	};


}

#endif /* REPEAT_VECS_H_ */
