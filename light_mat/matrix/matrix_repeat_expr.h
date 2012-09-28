/**
 * @file matrix_repeat_expr.h
 *
 * Expressions of repeating matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_REPEAT_EXPR_H_
#define LIGHTMAT_MATRIX_REPEAT_EXPR_H_

#include <light_mat/matrix/matrix_expr_base.h>
#include <light_mat/matrix/matrix_shape.h>

namespace lmat
{

	/********************************************
	 *
	 *  Expression classes
	 *
	 ********************************************/

	// traits

	template<typename Arg_HP, class Arg, int N>
	struct matrix_traits<horizontal_repeat_expr<Arg_HP, Arg, N> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value;
		static const int compile_time_num_cols = ct_cols<Arg>::value * N;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};

	template<typename Arg_HP, class Arg, int M>
	struct matrix_traits<vertical_repeat_expr<Arg_HP, Arg, M> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<Arg>::value * M;
		static const int compile_time_num_cols = ct_cols<Arg>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	// class definitions

	template<typename Arg_HP, class Arg, int N>
	class horizontal_repeat_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
		horizontal_repeat_expr<Arg_HP, Arg, N>,
		typename matrix_traits<Arg>::value_type>
	{
		typedef unary_expr_base<Arg_HP, Arg> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:

		LMAT_ENSURE_INLINE
		horizontal_repeat_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd, const index_t n)
		: base_t(arg_fwd), m_n(n)
		{
			check_arg( is_column(this->arg()), "repeat_col: the input argument should be a column." );
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
			return this->arg().nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_n.dim();
		}

	private:
		single_dim<N> m_n;
	};


	// vertical_repeat_expr

	template<typename Arg_HP, class Arg, int M>
	class vertical_repeat_expr
	: public unary_expr_base<Arg_HP, Arg>
	, public IMatrixXpr<
		vertical_repeat_expr<Arg_HP, Arg, M>,
		typename matrix_traits<Arg>::value_type>
	{
		typedef unary_expr_base<Arg_HP, Arg> base_t;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<Arg>::value, "Arg must be a matrix expression class.");
#endif

	public:

		LMAT_ENSURE_INLINE
		vertical_repeat_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd, const index_t m)
		: base_t(arg_fwd), m_m(m)
		{
			check_arg( is_row(this->arg()), "repeat_row: the input argument should be a row." );
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
			return m_m.dim();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return this->arg().ncolumns();
		}

	private:
		single_dim<M> m_m;
	};


	/********************************************
	 *
	 *  Expression mapping construction
	 *
	 ********************************************/

	// spec

	template<int N>
	struct hrep_t
	{
		LMAT_ENSURE_INLINE
		hrep_t() { }

		LMAT_ENSURE_INLINE
		hrep_t(const index_t n_)
		{
			check_arg( n_ == N, "The input arg to hrep_t should be equal to N" );
		}

		LMAT_ENSURE_INLINE
		index_t get_n() const { return N; }
	};

	template<>
	struct hrep_t<0>
	{
		const index_t n;

		LMAT_ENSURE_INLINE
		hrep_t(const index_t n_) : n(n_) { }

		LMAT_ENSURE_INLINE
		index_t get_n() const { return n; }
	};

	template<int M>
	struct vrep_t
	{
		LMAT_ENSURE_INLINE
		vrep_t() { }

		LMAT_ENSURE_INLINE
		vrep_t(const index_t m_)
		{
			check_arg( m_ == M, "The input arg to vrep_t should be equal to N" );
		}

		LMAT_ENSURE_INLINE
		index_t get_m() const { return M; }
	};

	template<>
	struct vrep_t<0>
	{
		const index_t m;

		LMAT_ENSURE_INLINE
		vrep_t(const index_t m_) : m(m_) { }

		LMAT_ENSURE_INLINE
		index_t get_m() const { return m; }
	};


	// maps

	template<class Arg, int N>
	struct unary_expr_verifier<hrep_t<N>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};

	template<class Arg, int M>
	struct unary_expr_verifier<vrep_t<M>, Arg>
	{
		static const bool value = is_mat_xpr<Arg>::value;
	};


	template<typename Arg_HP, class Arg, int N>
	struct unary_expr_map<hrep_t<N>, Arg_HP, Arg>
	{
		typedef horizontal_repeat_expr<Arg_HP, Arg, N> type;

		LMAT_ENSURE_INLINE
		static type get(const hrep_t<N>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd, spec.get_n());
		}
	};

	template<typename Arg_HP, class Arg, int M>
	struct unary_expr_map<vrep_t<M>, Arg_HP, Arg>
	{
		typedef vertical_repeat_expr<Arg_HP, Arg, M> type;

		LMAT_ENSURE_INLINE
		static type get(const vrep_t<M>& spec,
				const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		{
			return type(arg_fwd, spec.get_m());
		}
	};

	template<typename Arg_HP, typename T, int Mc, int Nc, int N>
	struct unary_expr_map<hrep_t<N>, Arg_HP, const_matrix<T, Mc, Nc> >
	{
		typedef const_matrix<T, Mc, N> type;

		LMAT_ENSURE_INLINE
		static type get(const hrep_t<N>& spec,
				const arg_forwarder<Arg_HP, const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(arg.nrows(), spec.get_n(), arg.value());
		}
	};


	template<typename Arg_HP, typename T, int Mc, int Nc, int M>
	struct unary_expr_map<vrep_t<M>, Arg_HP, const_matrix<T, Mc, Nc> >
	{
		typedef const_matrix<T, M, Nc> type;

		LMAT_ENSURE_INLINE
		static type get(const vrep_t<M>& spec,
				const arg_forwarder<Arg_HP, const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(spec.get_m(), arg.ncolumns(), arg.value());
		}
	};

	// construction functions

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<hrep_t<0>, ref_arg_t, Arg>::type
	hrep(const IMatrixXpr<Arg, T>& arg, const index_t n)
	{
		return make_expr(hrep_t<0>(n), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<hrep_t<N>, ref_arg_t, Arg>::type
	hrep(const IMatrixXpr<Arg, T>& arg, fix_int<N>)
	{
		return make_expr(hrep_t<N>(), ref_arg(arg.derived()));
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<vrep_t<0>, ref_arg_t, Arg>::type
	vrep(const IMatrixXpr<Arg, T>& arg, const index_t m)
	{
		return make_expr(vrep_t<0>(m), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int M>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<vrep_t<M>, ref_arg_t, Arg>::type
	vrep(const IMatrixXpr<Arg, T>& arg, fix_int<M>)
	{
		return make_expr(vrep_t<M>(), ref_arg(arg.derived()));
	}


}

#endif /* REPEAT_VECS_H_ */
