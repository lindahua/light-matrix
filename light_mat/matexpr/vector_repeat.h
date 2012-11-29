/**
 * @file vector_repeat.h
 *
 * Expression & evaluation of vector repeating
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_VECTOR_REPEAT_H_
#define LIGHTMAT_VECTOR_REPEAT_H_

#include <light_mat/matexpr/matrix_ewise_eval.h>
#include "bits/repeat_vecs_internal.h"

namespace lmat
{
	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	// expression traits

	template<typename Arg_HP, class Arg, int N>
	struct matrix_traits<repeat_col_expr<Arg_HP, Arg, N> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value;
		static const int ct_num_cols = ct_cols<Arg>::value * N;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Arg_HP, class Arg, int M>
	struct matrix_traits<repeat_row_expr<Arg_HP, Arg, M> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = ct_rows<Arg>::value * M;
		static const int ct_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	// expression

	template<typename Arg_HP, class Arg, int N>
	class repeat_col_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
		repeat_col_expr<Arg_HP, Arg, N>,
		typename matrix_traits<Arg>::value_type>
	{
		typedef unary_expr_base<Arg_HP, Arg> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:
		LMAT_ENSURE_INLINE
		repeat_col_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd, const index_t n)
		: base_t(arg_fwd), m_ncols(n)
		{
			check_arg( is_column(this->arg()),
					"repeat_col: the input argument should be a column." );
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * m_ncols.value();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols.value();
		}

	private:
		dimension<N> m_ncols;
	};


	template<typename Arg_HP, class Arg, int M>
	class repeat_row_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
		repeat_row_expr<Arg_HP, Arg, M>,
		typename matrix_traits<Arg>::value_type>
	{
		typedef unary_expr_base<Arg_HP, Arg> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:

		LMAT_ENSURE_INLINE
		repeat_row_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd, const index_t m)
		: base_t(arg_fwd), m_nrows(m)
		{
			check_arg( is_row(this->arg()), "repeat_row: the input argument should be a row." );
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * m_nrows.value();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_nrows.value();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		dimension<M> m_nrows;
	};


	/********************************************
	 *
	 *  Expression specification
	 *
	 ********************************************/

	// spec

	template<int N>
	struct repcol_t
	{
		LMAT_ENSURE_INLINE
		repcol_t() { }

		LMAT_ENSURE_INLINE
		repcol_t(const index_t n_)
		{
			check_arg( n_ == N, "The input arg to repcol_t should be equal to N" );
		}

		LMAT_ENSURE_INLINE
		index_t get_ncols() const { return N; }
	};

	template<>
	struct repcol_t<0>
	{
		const index_t ncols;

		LMAT_ENSURE_INLINE
		repcol_t() : ncols(0) { }

		LMAT_ENSURE_INLINE
		repcol_t(const index_t n_) : ncols(n_) { }

		LMAT_ENSURE_INLINE
		index_t get_ncols() const { return ncols; }
	};

	template<int M>
	struct reprow_t
	{
		LMAT_ENSURE_INLINE
		reprow_t() { }

		LMAT_ENSURE_INLINE
		reprow_t(const index_t m_)
		{
			check_arg( m_ == M, "The input arg to reprow_t should be equal to N" );
		}

		LMAT_ENSURE_INLINE
		index_t get_nrows() const { return M; }
	};

	template<>
	struct reprow_t<0>
	{
		const index_t nrows;

		LMAT_ENSURE_INLINE
		reprow_t() : nrows(0) { }

		LMAT_ENSURE_INLINE
		reprow_t(const index_t m_) : nrows(m_) { }

		LMAT_ENSURE_INLINE
		index_t get_nrows() const { return nrows; }
	};


	/********************************************
	 *
	 *  Expression verification and mapping
	 *
	 ********************************************/

	template<class Arg, int N>
	struct unary_expr_verifier<repcol_t<N>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<class Arg, int M>
	struct unary_expr_verifier<reprow_t<M>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};


	template<typename Arg_HP, class Arg, int N>
	struct unary_expr_map<repcol_t<N>, Arg_HP, Arg>
	{
		typedef repeat_col_expr<Arg_HP, Arg, N> type;

		LMAT_ENSURE_INLINE
		static type get(const repcol_t<N>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd, spec.get_ncols());
		}
	};

	template<typename Arg_HP, class Arg, int M>
	struct unary_expr_map<reprow_t<M>, Arg_HP, Arg>
	{
		typedef repeat_row_expr<Arg_HP, Arg, M> type;

		LMAT_ENSURE_INLINE
		static type get(const reprow_t<M>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd, spec.get_nrows());
		}
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<repcol_t<0>, ref_arg_t, Arg>::type
	rep_col(const IMatrixXpr<Arg, T>& arg, const index_t n)
	{
		return make_expr(repcol_t<0>(n), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<repcol_t<N>, ref_arg_t, Arg>::type
	rep_col(const IMatrixXpr<Arg, T>& arg, fix_int<N>)
	{
		return make_expr(repcol_t<N>(), ref_arg(arg.derived()));
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<reprow_t<0>, ref_arg_t, Arg>::type
	rep_row(const IMatrixXpr<Arg, T>& arg, const index_t m)
	{
		return make_expr(reprow_t<0>(m), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int M>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<reprow_t<M>, ref_arg_t, Arg>::type
	rep_row(const IMatrixXpr<Arg, T>& arg, fix_int<M>)
	{
		return make_expr(reprow_t<M>(), ref_arg(arg.derived()));
	}


	/********************************************
	 *
	 *  default evaluation
	 *
	 ********************************************/

	template<int M, int N>
	struct repcol_scheme
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		repcol_scheme(index_t m, index_t n)
		: shape(m, n) { }

		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const repeat_col_expr<Arg_HP, Arg, N>& sexpr, DMat& dmat)
		{
			internal::repcol_evaluate(shape, sexpr.arg(), dmat);
		}
	};

	template<int M, int N>
	struct reprow_scheme
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		reprow_scheme(index_t m, index_t n)
		: shape(m, n) { }

		template<typename Arg_HP, class Arg, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const repeat_row_expr<Arg_HP, Arg, M>& sexpr, DMat& dmat)
		{
			internal::reprow_evaluate(shape, sexpr.arg(), dmat);
		}
	};


	template<typename T, typename Arg_HP, class Arg, int N, class DMat>
	LMAT_ENSURE_INLINE
	inline repcol_scheme<common_ctrows<Arg, DMat>::value, N>
	get_default_eval_scheme(
			const repeat_col_expr<Arg_HP, Arg, N>& sexpr,
			IDenseMatrix<DMat, T>& dmat)
	{
		const int M = common_ctrows<Arg, DMat>::value;
		return repcol_scheme<M, N>(dmat.nrows(), dmat.ncolumns());
	}

	template<typename T, typename Arg_HP, class Arg, int M, class DMat>
	LMAT_ENSURE_INLINE
	inline reprow_scheme<M, common_ctcols<Arg, DMat>::value>
	get_default_eval_scheme(
			const repeat_row_expr<Arg_HP, Arg, M>& sexpr,
			IDenseMatrix<DMat, T>& dmat)
	{
		const int N = common_ctcols<Arg, DMat>::value;
		return reprow_scheme<M, N>(dmat.nrows(), dmat.ncolumns());
	}


	/********************************************
	 *
	 *  Percol accessors
	 *
	 ********************************************/

	// forward declaration

	template<typename Ker, class Arg>
	class repcol_accessor;

	template<typename Ker, typename Arg>
	class reprow_accessor;

	// class definitions

	template<class Arg>
	struct percol_macc_state_map<
		repcol_accessor<scalar_kernel_t, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef const arg_value_type* type;
	};

	template<class Arg>
	struct percol_macc_state_map<
		reprow_accessor<scalar_kernel_t, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef arg_value_type type;
	};


	template<class Arg>
	class repcol_accessor<scalar_kernel_t, Arg>
	: public IPerColMatrixScalarAccessor<
	  	  repcol_accessor<scalar_kernel_t, Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type value_type;

		template<typename Arg_HP, int N>
		repcol_accessor(const repeat_col_expr<Arg_HP, Arg, N>& expr)
		: m_cap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t i, const value_type* s) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE
		const value_type* col_state(const index_t j) const
		{
			return m_cap.pdata;
		}

	private:
		internal::repcol_cap<value_type, ct_rows<Arg>::value, is_dense_mat<Arg>::value> m_cap;

	};


	template<class Arg>
	class reprow_accessor<scalar_kernel_t, Arg>
	: public IPerColMatrixScalarAccessor<
	  	  reprow_accessor<scalar_kernel_t, Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type value_type;

		template<typename Arg_HP, int N>
		reprow_accessor(const repeat_row_expr<Arg_HP, Arg, N>& expr)
		: m_cap(expr.arg()) { }

		LMAT_ENSURE_INLINE
		value_type get_scalar(const index_t, const value_type& s) const
		{
			return s;
		}

		LMAT_ENSURE_INLINE
		value_type col_state(const index_t j) const
		{
			return m_cap[j];
		}

	private:
		internal::reprow_cap<Arg, ct_supports_linear_index_ex<Arg>::value> m_cap;

	};



	/********************************************
	 *
	 *  Accessor map & Cost model
	 *
	 ********************************************/

	template<typename Arg_HP, class Arg, int N, typename Ker>
	struct macc_accessor_map<repeat_col_expr<Arg_HP, Arg, N>, percol_macc, Ker>
	{
		typedef repcol_accessor<Ker, Arg> type;
	};

	template<typename Arg_HP, class Arg, int M, typename Ker>
	struct macc_accessor_map<repeat_row_expr<Arg_HP, Arg, M>, percol_macc, Ker>
	{
		typedef reprow_accessor<Ker, Arg> type;
	};


	template<typename Arg_HP, class Arg, int N, typename Ker>
	struct macc_cost<repeat_col_expr<Arg_HP, Arg, N>, linear_macc, Ker>
	{
		static const int value = MACC_CACHE_COST;
	};

	template<typename Arg_HP, class Arg, int M, typename Ker>
	struct macc_cost<repeat_row_expr<Arg_HP, Arg, M>, linear_macc, Ker>
	{
		static const int value = MACC_CACHE_COST;
	};

	template<typename Arg_HP, class Arg, int N, typename Ker>
	struct macc_cost<repeat_col_expr<Arg_HP, Arg, N>, percol_macc, Ker>
	{
		static const int value = 0;
	};

	template<typename Arg_HP, class Arg, int M, typename Ker>
	struct macc_cost<repeat_row_expr<Arg_HP, Arg, M>, percol_macc, Ker>
	{
		static const int value = 0;
	};

}


#endif


