/**
 * @file matrix_transpose_baseal.h
 *
 * Internal implementation of matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_INTERNAL_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_INTERNAL_H_

#include <light_mat/matrix/matrix_expr_base.h>

namespace lmat { namespace detail {

	// forward declaration

	template<class Col> class contcol_transpose_base;
	template<class Row> class controw_transpose_base;
	template<class Row> class regular_row_transpose_base;
	template<class Mat> class dense_transpose_base;
	template<class Expr> class generic_transpose_base;

	// traits

	template<class B>
	struct transpose_base_traits
	{
		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = false;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = false;
	};

	template<class Col>
	struct transpose_base_traits<contcol_transpose_base<Col> >
	{
		typedef typename unwrapped_expr<Col>::type arg_type;

		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = lmat::is_base_aligned<arg_type>::value;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = true;
	};

	template<class Row>
	struct transpose_base_traits<controw_transpose_base<Row> >
	{
		typedef typename unwrapped_expr<Row>::type arg_type;

		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = lmat::is_base_aligned<arg_type>::value;
		static const bool is_percol_aligned = is_base_aligned;
		static const bool is_linear_accessible = true;
	};

	template<class Row>
	struct transpose_base_traits<regular_row_transpose_base<Row> >
	{
		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = false;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = true;
	};

	/********************************************
	 *
	 *  transpose continuous vector
	 *
	 ********************************************/

	template<class Col>
	class contcol_transpose_base
	: public IDenseMatrix<contcol_transpose_base<Col>, typename matrix_traits<Col>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Col>::value_type )

		typedef typename unwrapped_expr<Col>::type arg_type;

	public:
		LMAT_ENSURE_INLINE contcol_transpose_base(const Col& col)
		: m_col(col)
		{
		}

		LMAT_ENSURE_INLINE const arg_type& arg() const
		{
			return m_col.get();
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t lead_dim() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return arg().ptr_data() + j;
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type, const index_type j) const
		{
			return ptr_data()[j];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return ptr_data()[idx];
		}

	private:
		obj_wrapper<Col> m_col;
	};


	template<class Row>
	class controw_transpose_base
	: public IDenseMatrix<controw_transpose_base<Row>, typename matrix_traits<Row>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Row>::value_type )

		typedef typename unwrapped_expr<Row>::type arg_type;

	public:
		LMAT_ENSURE_INLINE controw_transpose_base(const Row& row)
		: m_row(row)
		{
		}

		LMAT_ENSURE_INLINE const arg_type& arg() const
		{
			return m_row.get();
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_t lead_dim() const
		{
			return nrows();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type) const
		{
			return arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type) const
		{
			return ptr_data()[i];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return ptr_data()[idx];
		}

	private:
		obj_wrapper<Row> m_row;
	};


	/********************************************
	 *
	 *  transpose regular row
	 *
	 ********************************************/

	template<class Row>
	class regular_row_transpose_base
	: public IMatrixView<regular_row_transpose_base<Row>, typename matrix_traits<Row>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Row>::value_type )

		typedef typename unwrapped_expr<Row>::type arg_type;

	public:
		LMAT_ENSURE_INLINE regular_row_transpose_base(const Row& row)
		: m_row(row)
		{
		}

		LMAT_ENSURE_INLINE const arg_type& arg() const
		{
			return m_row.get();
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type) const
		{
			return arg().elem(0, i);
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return arg().elem(0, idx);
		}

	private:
		obj_wrapper<Row> m_row;
	};


	/********************************************
	 *
	 *  generic transpose
	 *
	 ********************************************/

	template<class Mat>
	class dense_transpose_base
	: public IMatrixXpr<dense_transpose_base<Mat>, typename matrix_traits<Mat>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Mat>::value_type )

		typedef typename unwrapped_expr<Mat>::type arg_type;

	public:
		LMAT_ENSURE_INLINE dense_transpose_base(const Mat& expr)
		: m_mat(expr)
		{
		}

		LMAT_ENSURE_INLINE const arg_type& arg() const
		{
			return m_mat.get();
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return arg().nrows();
		}

	private:
		obj_wrapper<Mat> m_mat;
	};


	template<class Expr>
	class generic_transpose_base
	: public IMatrixXpr<generic_transpose_base<Expr>, typename matrix_traits<Expr>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Expr>::value_type )

		typedef typename unwrapped_expr<Expr>::type arg_type;

	public:
		LMAT_ENSURE_INLINE generic_transpose_base(const Expr& expr)
		: m_expr(expr)
		{
		}

		LMAT_ENSURE_INLINE const arg_type& arg() const
		{
			return m_expr.get();
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return arg().nrows();
		}

	private:
		obj_wrapper<Expr> m_expr;
	};


	/********************************************
	 *
	 *  base type map
	 *
	 ********************************************/

	template<typename Expr_>
	struct matrix_transpose_base_map
	{
		typedef typename unwrapped_expr<Expr_>::type Expr;

		typedef typename
				if_<is_dense_mat<Expr>,
					typename
					if_<ct_is_col<Expr>,
						// is column
						contcol_transpose_base<Expr_>,
						typename
						if_<ct_is_row<Expr>,
							// is row
							typename
							if_<ct_has_continuous_layout<Expr>,
								controw_transpose_base<Expr_>,
								regular_row_transpose_base<Expr_>
							>::type,
							// is dense matrix (non-vector)
							dense_transpose_base<Expr_>
						>::type
					>::type,
					// generic expression
					generic_transpose_base<Expr_>
				>::type type;
	};


} }

#endif /* MATRIX_TRANSPOSE_INTERNAL_H_ */
