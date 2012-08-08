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

	template<typename Arg_HP, class Arg> class contcol_transpose_base;
	template<typename Arg_HP, class Arg> class controw_transpose_base;
	template<typename Arg_HP, class Arg> class regular_row_transpose_base;
	template<typename Arg_HP, class Arg> class dense_transpose_base;

	template<typename Arg_HP, class Arg> class colxpr_transpose_base;
	template<typename Arg_HP, class Arg> class rowxpr_transpose_base;
	template<typename Arg_HP, class Arg> class generic_transpose_base;

	// traits

	template<class B>
	struct transpose_base_traits
	{
		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = false;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = false;
	};

	template<typename Arg_HP, class Arg>
	struct transpose_base_traits<contcol_transpose_base<Arg_HP, Arg> >
	{
		static const bool has_continuous_layout = true;
		static const bool is_base_aligned = lmat::is_base_aligned<Arg>::value;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = true;
	};

	template<typename Arg_HP, class Arg>
	struct transpose_base_traits<controw_transpose_base<Arg_HP, Arg> >
	{
		static const bool has_continuous_layout = true;
		static const bool is_base_aligned = lmat::is_base_aligned<Arg>::value;
		static const bool is_percol_aligned = is_base_aligned;
		static const bool is_linear_accessible = true;
	};

	template<typename Arg_HP, class Arg>
	struct transpose_base_traits<regular_row_transpose_base<Arg_HP, Arg> >
	{
		static const bool has_continuous_layout = false;
		static const bool is_base_aligned = false;
		static const bool is_percol_aligned = false;
		static const bool is_linear_accessible = true;
	};


	// contcol_transpose_base

	template<typename Arg_HP, class Arg>
	class contcol_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IDenseMatrix<
	  	  transpose_expr<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE contcol_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t lead_dim() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return this->arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return this->arg().ptr_data() + j;
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type, const index_type j) const
		{
			return ptr_data()[j];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return ptr_data()[idx];
		}

	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const contcol_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		copy(s.ptr_data(), dst.derived());
	}


	// controw_transpose_base

	template<typename Arg_HP, class Arg>
	class controw_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IDenseMatrix<
	  	  transpose_expr<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE controw_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return this->arg().ncolumns();
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
			return this->arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type) const
		{
			return this->arg().ptr_data();
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type) const
		{
			return this->ptr_data()[i];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return this->ptr_data()[idx];
		}
	};


	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const controw_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		copy(s.ptr_data(), dst.derived());
	}


	// regular_row_transpose_base

	template<typename Arg_HP, class Arg>
	class regular_row_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixView<
	  	  transpose_expr<Arg_HP, Arg>,
	  	  typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE regular_row_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return this->arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type) const
		{
			return this->arg().elem(0, i);
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_t idx) const
		{
			return this->arg().elem(0, idx);
		}
	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const regular_row_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		const index_t m = s.nrows();
		typename matrix_traits<Arg>::value_type *pd = dst.ptr_data();
		const Arg& arg = s.arg();

		for (index_t i = 0; i < m; ++i)
		{
			pd[i] = arg.elem(0, i);
		}
	}


	// dense_transpose_base

	template<typename Arg_HP, class Arg>
	class dense_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<transpose_expr<Arg_HP, Arg>,
		typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE dense_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return this->arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return this->arg().nrows();
		}
	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const dense_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		const Arg& arg = s.arg();

		detail::transpose(arg.nrows(), arg.ncolumns(),
				arg.ptr_data(), arg.lead_dim(), dst.ptr_data(), dst.lead_dim());
	}


	// colxpr_transpose_base

	template<typename Arg_HP, class Arg>
	class colxpr_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<transpose_expr<Arg_HP, Arg>,
		typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE colxpr_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return 1;
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return this->arg().nrows();
		}
	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const colxpr_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		typedef typename matrix_traits<Arg>::value_type T;
		const int Len = binary_ctdim<ct_rows<Arg>::value, ct_cols<DMat>::value>::value;

		const Arg& arg = s.arg();

		if (has_continuous_layout(dst))
		{
			ref_matrix<T, Len, 1> dview(dst.ptr_data(), arg.nrows(), 1);
			default_evaluate(arg, dview);
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


	// rowxpr_transpose_base

	template<typename Arg_HP, class Arg>
	class rowxpr_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<transpose_expr<Arg_HP, Arg>,
		typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE rowxpr_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return this->arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return 1;
		}
	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const rowxpr_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		typedef typename matrix_traits<Arg>::value_type T;
		const int Len = binary_ctdim<ct_cols<Arg>::value, ct_rows<DMat>::value>::value;

		const Arg& arg = s.arg();
		ref_matrix<T, 1, Len> dview(dst.ptr_data(), 1, arg.ncolumns());
		default_evaluate(arg, dview);
	}



	/********************************************
	 *
	 *  generic transpose
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg>
	class generic_transpose_base
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<transpose_expr<Arg_HP, Arg>,
		typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef unary_expr_base<Arg_HP, Arg> base_t;
		LMAT_MAT_TRAITS_CDEFS( typename matrix_traits<Arg>::value_type )

	public:
		LMAT_ENSURE_INLINE generic_transpose_base(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}

		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return this->arg().nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return this->arg().size();
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return this->arg().ncolumns();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return this->arg().nrows();
		}
	};

	template<typename Arg_HP, class Arg, class DMat>
	LMAT_ENSURE_INLINE
	inline void transbase_evaluate_to(const generic_transpose_base<Arg_HP, Arg>& s,
			IDenseMatrix<DMat, typename matrix_traits<Arg>::value_type>& dst)
	{
		dense_matrix<typename matrix_traits<Arg>::value_type,
			ct_rows<Arg>::value,
			ct_cols<Arg>::value> mat = s.arg();

		detail::transpose(mat.nrows(), mat.ncolumns(),
				mat.ptr_data(), mat.lead_dim(), dst.ptr_data(), dst.lead_dim());
	}


	/********************************************
	 *
	 *  base type map
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg>
	struct matrix_transpose_base_map
	{
		typedef typename
				if_<is_dense_mat<Arg>,
					// is dense
					typename
					if_<ct_is_col<Arg>,
						// is column
						contcol_transpose_base<Arg_HP, Arg>,
						typename
						if_<ct_is_row<Arg>,
							// is row
							typename
							if_<ct_has_continuous_layout<Arg>,
								controw_transpose_base<Arg_HP, Arg>,
								regular_row_transpose_base<Arg_HP, Arg>
							>::type,
							// is dense matrix (non-vector)
							dense_transpose_base<Arg_HP, Arg>
						>::type
					>::type,
					// non dense
					typename
					if_<ct_is_row<Arg>,
						// is_row
						rowxpr_transpose_base<Arg_HP, Arg>,
						typename
						if_<ct_is_col<Arg>,
							// is column
							colxpr_transpose_base<Arg_HP, Arg>,
							// generic non-vector
							generic_transpose_base<Arg_HP, Arg>
						>::type
					>::type
				>::type type;
	};


	/********************************************
	 *
	 *  reflection
	 *
	 ********************************************/

#define LMAT_DEFINE_TRANS_BASE_TYPENAME( tyname ) \
	template<typename Arg_HP, class Arg> \
	LMAT_ENSURE_INLINE \
	inline const char* trans_base_typename(const tyname##_transpose_base<Arg_HP, Arg>& ) \
	{ return #tyname; }

	LMAT_DEFINE_TRANS_BASE_TYPENAME( contcol )
	LMAT_DEFINE_TRANS_BASE_TYPENAME( controw )
	LMAT_DEFINE_TRANS_BASE_TYPENAME( regular_row )
	LMAT_DEFINE_TRANS_BASE_TYPENAME( dense )

	LMAT_DEFINE_TRANS_BASE_TYPENAME( colxpr )
	LMAT_DEFINE_TRANS_BASE_TYPENAME( rowxpr )
	LMAT_DEFINE_TRANS_BASE_TYPENAME( generic )


	/********************************************
	 *
	 *  transpose_expr classes
	 *
	 ********************************************/

	// traits

	template<typename Arg_HP, class Arg>
	struct matrix_traits<transpose_expr<Arg_HP, Arg> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_cols<Arg>::value;
		static const int compile_time_num_cols = ct_rows<Arg>::value;

		static const bool is_readonly = true;
		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Arg_HP, class Arg>
	struct ct_has_continuous_layout<transpose_expr<Arg_HP, Arg> >
	{
		typedef typename matrix_transpose_base_map<Arg_HP, Arg>::type base_t;
		static const bool value = transpose_base_traits<base_t>::has_continuous_layout;
	};

	template<typename Arg_HP, class Arg>
	struct is_base_aligned<transpose_expr<Arg_HP, Arg> >
	{
		typedef typename matrix_transpose_base_map<Arg_HP, Arg>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_base_aligned;
	};

	template<typename Arg_HP, class Arg>
	struct is_percol_aligned<transpose_expr<Arg_HP, Arg> >
	{
		typedef typename matrix_transpose_base_map<Arg_HP, Arg>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_percol_aligned;
	};

	template<typename Arg_HP, class Arg>
	struct is_linear_accessible<transpose_expr<Arg_HP, Arg> >
	{
		typedef typename matrix_transpose_base_map<Arg_HP, Arg>::type base_t;
		static const bool value = transpose_base_traits<base_t>::is_linear_accessible;
	};


	// class

	template<typename Arg_HP, class Arg>
	class transpose_expr : public matrix_transpose_base_map<Arg_HP, Arg>::type
	{
		typedef typename matrix_transpose_base_map<Arg_HP, Arg>::type base_t;

	public:
		LMAT_ENSURE_INLINE
		transpose_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: base_t(arg_fwd)
		{
		}
	};


	template<class Arg>
	struct unary_expr_verifier<transpose_t, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<typename Arg_HP, class Arg>
	struct unary_expr_map<transpose_t, Arg_HP, Arg>
	{
		typedef transpose_expr<Arg_HP, Arg> type;

		LMAT_ENSURE_INLINE
		static type get(transpose_t, const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd);
		}
	};


	// evaluation

	struct matrix_transpose_policy { };

	template<typename Arg_HP, class Arg, class Dst>
	struct default_matrix_eval_policy<transpose_expr<Arg_HP, Arg>, Dst>
	{
		typedef matrix_transpose_policy type;
	};

	template<typename Arg_HP, class Arg, class Dst>
	LMAT_ENSURE_INLINE
	inline void evaluate(
			const transpose_expr<Arg_HP, Arg>& expr,
			IDenseMatrix<Dst, typename matrix_traits<Arg>::value_type>& dst,
			matrix_transpose_policy)
	{
		transbase_evaluate_to(expr, dst.derived());
	}


}

#endif /* MATRIX_TRANSPOSE_H_ */
