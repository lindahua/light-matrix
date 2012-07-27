/**
 * @file matrix_transpose.h
 *
 * Matrix transposition
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_TRANSPOSE_H_
#define LIGHTMAT_MATRIX_TRANSPOSE_H_

#include <light_mat/matrix/matrix_expr_base.h>
#include "bits/mat_transpose_impl.h"

namespace lmat
{

	/********************************************
	 *
	 *  transpose_base classes
	 *
	 ********************************************/

	// forward declaration

	template<class Col> class contcol_transpose_base;
	template<class Row> class controw_transpose_base;
	template<class Row> class regular_row_transpose_base;
	template<class Mat> class dense_transpose_base;

	template<class Expr> class colxpr_transpose_base;
	template<class Expr> class rowxpr_transpose_base;
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


	template<class Col>
	class contcol_transpose_base
	: public IDenseMatrix<transpose_expr<Col>, typename matrix_traits<Col>::value_type>
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


	template<class Col, class DMat>
	LMAT_ENSURE_INLINE
	void base_evaluate_to(const contcol_transpose_base<Col>& s,
			IDenseMatrix<DMat, typename matrix_traits<Col>::value_type>& dst)
	{
		copy(s.ptr_data(), dst.derived());
	}



	template<class Row>
	class controw_transpose_base
	: public IDenseMatrix<transpose_expr<Row>, typename matrix_traits<Row>::value_type>
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


	template<class Row, class DMat>
	LMAT_ENSURE_INLINE
	void base_evaluate_to(const controw_transpose_base<Row>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		copy(s.ptr_data(), dst.derived());
	}



	/********************************************
	 *
	 *  transpose regular row
	 *
	 ********************************************/

	template<class Row>
	class regular_row_transpose_base
	: public IMatrixView<transpose_expr<Row>, typename matrix_traits<Row>::value_type>
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

	template<class Row, class DMat>
	LMAT_ENSURE_INLINE
	void base_evaluate_to(const regular_row_transpose_base<Row>& s,
			IDenseMatrix<DMat, typename matrix_traits<Row>::value_type>& dst)
	{
		const index_t m = s.nrows();
		typename matrix_traits<Row>::value_type *pd = dst.ptr_data();
		const typename unwrapped_expr<Row>::type& row = s.arg();

		for (index_t i = 0; i < m; ++i)
		{
			pd[i] = row.elem(0, i);
		}
	}



	/********************************************
	 *
	 *  dense transpose
	 *
	 ********************************************/

	template<class Mat>
	class dense_transpose_base
	: public IMatrixXpr<transpose_expr<Mat>, typename matrix_traits<Mat>::value_type>
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

	template<class Mat, class DMat>
	LMAT_ENSURE_INLINE
	void base_evaluate_to(const dense_transpose_base<Mat>& s,
			IDenseMatrix<DMat, typename matrix_traits<Mat>::value_type>& dst)
	{
		const typename unwrapped_expr<Mat>::type& mat = s.arg();

		detail::transpose(mat.nrows(), mat.ncolumns(),
				mat.ptr_data(), mat.lead_dim(), dst.ptr_data(), dst.lead_dim());
	}


	/********************************************
	 *
	 *  col expression transpose
	 *
	 ********************************************/

	template<class Expr>
	class colxpr_transpose_base
	: public IMatrixXpr<transpose_expr<Expr>, typename matrix_traits<Expr>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Expr>::value_type )

		typedef typename unwrapped_expr<Expr>::type arg_type;

	public:
		LMAT_ENSURE_INLINE colxpr_transpose_base(const Expr& expr)
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
			return 1;
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return arg().nrows();
		}

	private:
		obj_wrapper<Expr> m_expr;
	};

	template<class Expr, class DMat>
	inline
	void base_evaluate_to(const colxpr_transpose_base<Expr>& s,
			IDenseMatrix<DMat, typename matrix_traits<Expr>::value_type>& dst)
	{
		typedef typename matrix_traits<Expr>::value_type T;
		const int Len = binary_ctdim<ct_rows<Expr>::value, ct_cols<DMat>::value>::value;

		const typename unwrapped_expr<Expr>::type& arg = s.arg();

		if (has_continuous_layout(dst))
		{
			ref_matrix<T, Len, 1> dview(dst.ptr_data(), arg.nrows(), 1);
			evaluate_to(arg, dview);
		}
		else
		{
			dense_matrix<T, Len, 1> tmp(arg);
			const index_t n = arg.nrows();
			for (index_t i = 0; i < n; ++i)
			{
				dst.elem(0, i) = tmp[i];
			}
		}
	}


	/********************************************
	 *
	 *  row expression transpose
	 *
	 ********************************************/

	template<class Expr>
	class rowxpr_transpose_base
	: public IMatrixXpr<transpose_expr<Expr>, typename matrix_traits<Expr>::value_type>
	{
	public:
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Expr>::value_type )

		typedef typename unwrapped_expr<Expr>::type arg_type;

	public:
		LMAT_ENSURE_INLINE rowxpr_transpose_base(const Expr& expr)
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
			return 1;
		}

	private:
		obj_wrapper<Expr> m_expr;
	};

	template<class Expr, class DMat>
	inline
	void base_evaluate_to(const rowxpr_transpose_base<Expr>& s,
			IDenseMatrix<DMat, typename matrix_traits<Expr>::value_type>& dst)
	{
		typedef typename matrix_traits<Expr>::value_type T;
		const int Len = binary_ctdim<ct_cols<Expr>::value, ct_rows<DMat>::value>::value;

		const typename unwrapped_expr<Expr>::type& arg = s.arg();
		ref_matrix<T, 1, Len> dview(dst.ptr_data(), 1, arg.ncolumns());
		evaluate_to(arg, dview);
	}



	/********************************************
	 *
	 *  generic transpose
	 *
	 ********************************************/

	template<class Expr>
	class generic_transpose_base
	: public IMatrixXpr<transpose_expr<Expr>, typename matrix_traits<Expr>::value_type>
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

	template<class Expr, class DMat>
	LMAT_ENSURE_INLINE
	void base_evaluate_to(const generic_transpose_base<Expr>& s,
			IDenseMatrix<DMat, typename matrix_traits<Expr>::value_type>& dst)
	{
		dense_matrix<typename matrix_traits<Expr>::value_type,
			ct_rows<Expr>::value,
			ct_cols<Expr>::value> mat = s.arg();

		detail::transpose(mat.nrows(), mat.ncolumns(),
				mat.ptr_data(), mat.lead_dim(), dst.ptr_data(), dst.lead_dim());
	}


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
					// is dense
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
					// non dense
					typename
					if_<ct_is_row<Expr>,
						// is_row
						rowxpr_transpose_base<Expr_>,
						typename
						if_<ct_is_col<Expr>,
							// is column
							colxpr_transpose_base<Expr_>,
							// generic non-vector
							generic_transpose_base<Expr_>
						>::type
					>::type
				>::type type;
	};


	/********************************************
	 *
	 *  retrospection
	 *
	 ********************************************/

	template<class Col>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const contcol_transpose_base<Col>& )
	{
		return "contcol";
	}

	template<class Row>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const controw_transpose_base<Row>& )
	{
		return "controw";
	}

	template<class Row>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const regular_row_transpose_base<Row>& )
	{
		return "regular_row";
	}

	template<class Mat>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const dense_transpose_base<Mat>& )
	{
		return "dense";
	}

	template<class Expr>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const colxpr_transpose_base<Expr>& )
	{
		return "colxpr";
	}

	template<class Expr>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const rowxpr_transpose_base<Expr>& )
	{
		return "rowxpr";
	}

	template<class Expr>
	LMAT_ENSURE_INLINE
	const char* trans_base_typename(const generic_transpose_base<Expr>& )
	{
		return "generic";
	}


	/********************************************
	 *
	 *  transpose_expr classes
	 *
	 ********************************************/

	// traits

	template<class Expr>
	struct matrix_traits<transpose_expr<Expr> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_cols<Expr>::value;
		static const int compile_time_num_cols = ct_rows<Expr>::value;

		static const bool is_readonly = true;
		typedef typename matrix_traits<Expr>::value_type value_type;
	};

	template<class Expr>
	struct ct_has_continuous_layout<transpose_expr<Expr> >
	{
		typedef typename matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = transpose_base_traits<base_t>::has_continuous_layout;
	};

	template<class Expr>
	struct is_base_aligned<transpose_expr<Expr> >
	{
		typedef typename matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_base_aligned;
	};

	template<class Expr>
	struct is_percol_aligned<transpose_expr<Expr> >
	{
		typedef typename matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_percol_aligned;
	};

	template<class Expr>
	struct is_linear_accessible<transpose_expr<Expr> >
	{
		typedef typename matrix_transpose_base_map<Expr>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_linear_accessible;
	};

	// class

	template<class Expr>
	class transpose_expr : public matrix_transpose_base_map<Expr>::type
	{
		typedef typename matrix_transpose_base_map<Expr>::type base_t;

	public:
		LMAT_ENSURE_INLINE
		transpose_expr(const Expr& expr)
		: base_t(expr)
		{
		}
	};

	// evaluation

	template<class Expr, class DMat>
	LMAT_ENSURE_INLINE
	void evaluate_to(const transpose_expr<Expr>& s,
			IDenseMatrix<DMat, typename matrix_traits<Expr>::value_type>& dst)
	{
		base_evaluate_to(s, dst);
	}

}

#endif /* MATRIX_TRANSPOSE_H_ */
