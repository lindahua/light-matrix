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

#include <light_mat/matrix/matrix_classes.h>
#include "bits/repeat_vecs_internal.h"

namespace lmat
{
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

	public:
		typedef typename matrix_traits<Col>::value_type value_type;

		LMAT_ENSURE_INLINE
		repeat_col_expr(const Col& col, const index_t n)
		: m_col(col), m_ncols(n)
		{
			check_arg( is_column(col), "repeat_col: the input col is NOT a column." );
		}

		LMAT_ENSURE_INLINE const Col& column() const
		{
			return m_col;
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
			return m_col.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols.dim();
		}

	private:
		const Col& m_col;
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

	public:
		typedef typename matrix_traits<Row>::value_type value_type;

		LMAT_ENSURE_INLINE
		repeat_row_expr(const Row& row, const index_t m)
		: m_row(row), m_nrows(m)
		{
			check_arg( is_row(row), "repeat_row: the input row is NOT a row." );
		}

		LMAT_ENSURE_INLINE const Row& row() const
		{
			return m_row;
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
			return m_row.ncolumns();
		}

	private:
		const Row& m_row;
		single_dim<M> m_nrows;
	};


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
	inline typename repcol_type_map<const_matrix<T, Mc, Nc>, DynamicDim>::type
	repcol(const const_matrix<T, Mc, Nc>& col, const index_t n)
	{
		return const_matrix<T, Mc, DynamicDim>(col.nrows(), n, col.value());
	}

	template<typename T, int Mc, int Nc, int N>
	LMAT_ENSURE_INLINE
	inline typename repcol_type_map<const_matrix<T, Mc, Nc>, N>::type
	repcol(const const_matrix<T, Mc, Nc>& col, fixed_dim<N>)
	{
		return const_matrix<T, Mc, N>(col.nrows(), N, col.value());
	}

	template<typename T, int Mc, int Nc>
	LMAT_ENSURE_INLINE
	inline typename reprow_type_map<const_matrix<T, Mc, Nc>, DynamicDim>::type
	reprow(const const_matrix<T, Mc, Nc>& row, const index_t m)
	{
		return const_matrix<T, DynamicDim, Nc>(m, row.ncolumns(), row.value());
	}

	template<typename T, int Mc, int Nc, int M>
	LMAT_ENSURE_INLINE
	inline typename reprow_type_map<const_matrix<T, Mc, Nc>, M>::type
	reprow(const const_matrix<T, Mc, Nc>& row, fixed_dim<M>)
	{
		return const_matrix<T, M, Nc>(M, row.ncolumns(), row.value());
	}


	/********************************************
	 *
	 *  Expression evaluation
	 *
	 ********************************************/

	template<class Col, class DMat>
	inline void evaluate_to(const repeat_col_expr<Col, 1>& s,
			IDenseMatrix<DMat, typename matrix_traits<Col>::value_type>& dst)
	{
		evaluate_to(s.column(), dst);
	}


	template<class Col, int N, class DMat>
	inline void evaluate_to(const repeat_col_expr<Col, N>& s,
			IDenseMatrix<DMat, typename matrix_traits<Col>::value_type>& dst)
	{
		typedef typename matrix_traits<Col>::value_type T;

		if ( is_column(s) )
		{
			const int M = binary_ct_rows<Col, DMat>::value;
			ref_col<T, M> dview(dst.ptr_data(), dst.nrows());
			evaluate_to(s.column(), dview);
		}
		else
		{
			typedef typename detail::repcol_ewrapper_map<Col>::type wrapper_t;
			wrapper_t col_wrap(s.column());

			const index_t m = col_wrap.nrows();
			if (m == 1)
			{
				fill(dst, col_wrap[0]);
			}
			else
			{
				const index_t n = s.ncolumns();
				for (index_t j = 0; j < n; ++j)
				{
					copy_mem(m, col_wrap.data(), dst.ptr_col(j));
				}
			}
		}
	}

	template<class Row, class DMat>
	inline void evaluate_to(const repeat_row_expr<Row, 1>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		evaluate_to( s.row(), dst );
	}

	template<class Row, int M, class DMat>
	inline void evaluate_to(const repeat_row_expr<Row, M>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		typedef typename matrix_traits<Row>::value_type T;

		if ( is_row(s) )
		{
			const int N = binary_ct_cols<Row, DMat>::value;
			if (has_continuous_layout<DMat>::value)
			{
				ref_row<T, N> dview(dst.ptr_data(), dst.ncolumns());
				evaluate_to(s.row(), dview);
			}
			else
			{
				ref_matrix_ex<T, 1, N> dview(dst.ptr_data(), 1, dst.ncolumns(), dst.lead_dim());
				evaluate_to(s.row(), dview);
			}
		}
		else
		{
			typedef typename detail::reprow_ewrapper_map<Row>::type wrapper_t;
			wrapper_t row_wrap(s.row());

			const index_t n = row_wrap.ncolumns();

			if (M == 0)
			{
				const index_t m = s.nrows();
				if (n == 1)
				{
					fill_mem(m, dst.ptr_data(), row_wrap[0]);
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						fill_mem(m, dst.ptr_col(j), row_wrap[j]);
				}
			}
			else
			{
				if (n == 1)
				{
					fill_mem(M, dst.ptr_data(), row_wrap[0]);
				}
				else
				{
					for (index_t j = 0; j < n; ++j)
						fill_mem(M, dst.ptr_col(j), row_wrap[j]);
				}
			}
		}
	}

}

#endif /* REPEAT_VECS_H_ */
