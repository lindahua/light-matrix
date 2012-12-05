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

	template<class QArg, int N>
	struct matrix_traits<repeat_col_expr<QArg, N> >
	{
		typedef typename QArg::argument arg_t;

		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<arg_t>::value;
		static const int ct_num_cols = meta::ncols<arg_t>::value * N;

		static const bool is_readonly = true;

		typedef typename matrix_traits<arg_t>::value_type value_type;
		typedef typename matrix_traits<arg_t>::domain domain;
	};

	template<class QArg, int M>
	struct matrix_traits<repeat_row_expr<QArg, M> >
	{
		typedef typename QArg::argument arg_t;

		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<arg_t>::value * M;
		static const int ct_num_cols = meta::ncols<arg_t>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<arg_t>::value_type value_type;
		typedef typename matrix_traits<arg_t>::domain domain;
	};

	// expression

	template<class QArg, int N>
	class repeat_col_expr
	: public IMatrixXpr<
		repeat_col_expr<QArg, N>,
		typename matrix_traits<typename QArg::argument>::value_type>
	{
	public:
		LMAT_ENSURE_INLINE
		repeat_col_expr(const typename QArg::forwarder& arg_fwd, const index_t n)
		: arg_holder(arg_fwd), m_ncols(n)
		{
			check_arg( is_column(this->arg()),
					"repeat_col: the input argument should be a column." );
		}

		LMAT_ENSURE_INLINE const typename QArg::argument& arg() const
		{
			return arg_holder.get();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return nrows() * m_ncols.value();
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_ncols.value();
		}

	private:
		typename QArg::holder arg_holder;
		dimension<N> m_ncols;
	};


	template<class QArg, int M>
	class repeat_row_expr
	: public IMatrixXpr<
		repeat_row_expr<QArg, M>,
		typename matrix_traits<typename QArg::argument>::value_type>
	{
	public:

		LMAT_ENSURE_INLINE
		repeat_row_expr(const typename QArg::forwarder& arg_fwd, const index_t m)
		: arg_holder(arg_fwd), m_nrows(m)
		{
			check_arg( is_row(this->arg()), "repeat_row: the input argument should be a row." );
		}

		LMAT_ENSURE_INLINE const typename QArg::argument& arg() const
		{
			return arg_holder.get();
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
			return arg().ncolumns();
		}

	private:
		typename QArg::holder arg_holder;
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
	struct expr_verifier<repcol_t<N>, LMAT_TYPELIST_1(Arg) >
	{
		static const bool value = meta::is_mat_xpr<Arg>::value;
	};

	template<class Arg, int M>
	struct expr_verifier<reprow_t<M>, LMAT_TYPELIST_1(Arg) >
	{
		static const bool value = meta::is_mat_xpr<Arg>::value;
	};


	template<class QArg, int N>
	struct expr_map<repcol_t<N>, LMAT_TYPELIST_1(QArg) >
	{
		typedef repeat_col_expr<QArg, N> type;

		LMAT_ENSURE_INLINE
		static type get(const repcol_t<N>& spec,
				const tied_forwarder< LMAT_TYPELIST_1(QArg) >& tfwd)
		{
			return type(tfwd.arg1_fwd, spec.get_ncols());
		}
	};


	template<class QArg, int M>
	struct expr_map<reprow_t<M>, LMAT_TYPELIST_1(QArg) >
	{
		typedef repeat_row_expr<QArg, M> type;

		LMAT_ENSURE_INLINE
		static type get(const reprow_t<M>& spec,
				const tied_forwarder< LMAT_TYPELIST_1(QArg) > & tfwd)
		{
			return type(tfwd.arg1_fwd, spec.get_nrows());
		}
	};



	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename expr_map<repcol_t<0>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type
	rep_col(const IMatrixXpr<Arg, T>& arg, const index_t n)
	{
		typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;
		return make_expr(repcol_t<0>(n), tfwd_t(arg.derived()));
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	inline typename expr_map<repcol_t<N>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type
	rep_col(const IMatrixXpr<Arg, T>& arg, fix_int<N>)
	{
		typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;
		return make_expr(repcol_t<N>(), tfwd_t(arg.derived()));
	}


	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename expr_map<reprow_t<0>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type
	rep_row(const IMatrixXpr<Arg, T>& arg, const index_t n)
	{
		typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;
		return make_expr(reprow_t<0>(n), tfwd_t(arg.derived()));
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	inline typename expr_map<reprow_t<N>, LMAT_TYPELIST_1(CRefArg<Arg>) >::type
	rep_row(const IMatrixXpr<Arg, T>& arg, fix_int<N>)
	{
		typedef tied_forwarder< LMAT_TYPELIST_1(CRefArg<Arg>) > tfwd_t;
		return make_expr(reprow_t<N>(), tfwd_t(arg.derived()));
	}


	/********************************************
	 *
	 *  default evaluation
	 *
	 ********************************************/

	template<class QArg, int N, class DMat>
	LMAT_ENSURE_INLINE
	void evaluate(const repeat_col_expr<QArg, N>& sexpr,
			IDenseMatrix<DMat, typename matrix_traits<typename QArg::argument>::value_type>& dmat)
	{
		matrix_shape< meta::nrows<typename QArg::argument>::value, N> shape(
				common_nrows(sexpr, dmat),
				common_ncols(sexpr, dmat));

		internal::repcol_evaluate(shape, sexpr.arg(), dmat.derived());
	}

	template<class QArg, int M, class DMat>
	LMAT_ENSURE_INLINE
	void evaluate(const repeat_row_expr<QArg, M>& sexpr,
			IDenseMatrix<DMat, typename matrix_traits<typename QArg::argument>::value_type>& dmat)
	{
		matrix_shape<M, meta::ncols<typename QArg::argument>::value> shape(
				common_nrows(sexpr, dmat),
				common_ncols(sexpr, dmat));

		internal::reprow_evaluate(shape, sexpr.arg(), dmat.derived());
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
		repcol_accessor<scalar_ker, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef const arg_value_type* type;
	};

	template<class Arg>
	struct percol_macc_state_map<
		reprow_accessor<scalar_ker, Arg> >
	{
		typedef typename matrix_traits<Arg>::value_type arg_value_type;
		typedef arg_value_type type;
	};


	template<class Arg>
	class repcol_accessor<scalar_ker, Arg>
	: public IPerColMatrixScalarAccessor<
	  	  repcol_accessor<scalar_ker, Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type value_type;

		template<class QArg, int N>
		repcol_accessor(const repeat_col_expr<QArg, N>& expr)
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
		internal::repcol_cap<value_type, meta::nrows<Arg>::value, meta::is_dense_mat<Arg>::value> m_cap;

	};


	template<class Arg>
	class reprow_accessor<scalar_ker, Arg>
	: public IPerColMatrixScalarAccessor<
	  	  reprow_accessor<scalar_ker, Arg>, typename matrix_traits<Arg>::value_type>
	{
	public:
		typedef typename matrix_traits<Arg>::value_type value_type;

		template<class QArg, int N>
		reprow_accessor(const repeat_row_expr<QArg, N>& expr)
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
		internal::reprow_cap<Arg, meta::supports_linear_index_ex<Arg>::value> m_cap;

	};



	/********************************************
	 *
	 *  Accessor map & Cost model
	 *
	 ********************************************/

	template<class QArg, int N, typename Ker>
	struct macc_accessor_map<repeat_col_expr<QArg, N>, macc_policy<percol_macc, Ker> >
	{
		typedef repcol_accessor<Ker, typename QArg::argument> type;
	};

	template<class QArg, int M, typename Ker>
	struct macc_accessor_map<repeat_row_expr<QArg, M>, macc_policy<percol_macc, Ker> >
	{
		typedef reprow_accessor<Ker, typename QArg::argument> type;
	};


	template<class QArg, int N, typename Ker>
	struct macc_cost<repeat_col_expr<QArg, N>, linear_macc, Ker>
	{
		static const int value = MACC_CACHE_COST;
	};

	template<class QArg, int M, typename Ker>
	struct macc_cost<repeat_row_expr<QArg, M>, linear_macc, Ker>
	{
		static const int value = MACC_CACHE_COST;
	};

	template<class QArg, int N, typename Ker>
	struct macc_cost<repeat_col_expr<QArg, N>, percol_macc, Ker>
	{
		static const int value = 0;
	};

	template<class QArg, int M, typename Ker>
	struct macc_cost<repeat_row_expr<QArg, M>, percol_macc, Ker>
	{
		static const int value = 0;
	};

}


#endif


