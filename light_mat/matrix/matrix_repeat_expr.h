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

	template<class Arg_Holder, int N>
	struct matrix_traits<horizontal_repeat_expr<Arg_Holder, N> >
	{
		typedef typename Arg_Holder::arg_type arg_type;

		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<arg_type>::value;
		static const int compile_time_num_cols = ct_cols<arg_type>::value * N;

		static const bool is_readonly = true;

		typedef typename matrix_traits<arg_type>::value_type value_type;
		typedef typename matrix_traits<arg_type>::domain domain;
	};

	template<class Arg_Holder, int M>
	struct matrix_traits<vertical_repeat_expr<Arg_Holder, M> >
	{
		typedef typename Arg_Holder::arg_type arg_type;

		static const int num_dimensions = 2;
		static const int compile_time_num_rows = ct_rows<arg_type>::value * M;
		static const int compile_time_num_cols = ct_cols<arg_type>::value;

		static const bool is_readonly = true;

		typedef typename matrix_traits<arg_type>::value_type value_type;
		typedef typename matrix_traits<arg_type>::domain domain;
	};


	// class definitions

	template<class Arg_Holder, int N>
	class horizontal_repeat_expr
	: public unary_expr<Arg_Holder>
	, public IMatrixXpr<
		horizontal_repeat_expr<Arg_Holder, N>,
		typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
		typedef unary_expr<Arg_Holder> base_t;
		typedef typename base_t::arg_type arg_type;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<arg_type>::value, "Arg must be a matrix expression class.");
#endif

	public:
		LMAT_ENSURE_INLINE
		horizontal_repeat_expr(const typename base_t::arg_forwarder& arg_fwd, const index_t n)
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

	template<class Arg_Holder, int M>
	class vertical_repeat_expr
	: public unary_expr<Arg_Holder>
	, public IMatrixXpr<
		vertical_repeat_expr<Arg_Holder, M>,
		typename matrix_traits<typename Arg_Holder::arg_type>::value_type>
	{
		typedef unary_expr<Arg_Holder> base_t;
		typedef typename base_t::arg_type arg_type;

#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_mat_xpr<arg_type>::value, "Arg must be a matrix expression class.");
#endif

	public:
		LMAT_ENSURE_INLINE
		vertical_repeat_expr(const typename base_t::arg_forwarder& arg_fwd, const index_t m)
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
	struct hrep_t<DynamicDim>
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
	struct vrep_t<DynamicDim>
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


	template<class Arg_Holder, int N>
	struct unary_expr_map<hrep_t<N>, Arg_Holder>
	{
		typedef horizontal_repeat_expr<Arg_Holder, N> type;

		LMAT_ENSURE_INLINE
		static type get(const hrep_t<N>& spec,
				const typename arg_forwarder<Arg_Holder>::type& arg_fwd)
		{
			return type(arg_fwd, spec.get_n());
		}
	};

	template<class Arg_Holder, int M>
	struct unary_expr_map<vrep_t<M>, Arg_Holder>
	{
		typedef vertical_repeat_expr<Arg_Holder, M> type;

		LMAT_ENSURE_INLINE
		static type get(const vrep_t<M>& spec,
				const typename arg_forwarder<Arg_Holder>::type& arg_fwd)
		{
			return type(arg_fwd, spec.get_m());
		}
	};

	template<typename T, int Mc, int Nc, int N>
	struct unary_expr_map<hrep_t<N>, ref_arg_holder<const_matrix<T, Mc, Nc> > >
	{
		typedef const_matrix<T, Mc, N> type;

		LMAT_ENSURE_INLINE
		static type get(const hrep_t<N>& spec,
				const ref_arg_forwarder<const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(arg.nrows(), spec.get_n(), arg.value());
		}
	};

	template<typename T, int Mc, int Nc, int N>
	struct unary_expr_map<hrep_t<N>, copy_arg_holder<const_matrix<T, Mc, Nc> > >
	{
		typedef const_matrix<T, Mc, N> type;

		LMAT_ENSURE_INLINE
		static type get(const hrep_t<N>& spec,
				const copy_arg_forwarder<const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(arg.nrows(), spec.get_n(), arg.value());
		}
	};

	template<typename T, int Mc, int Nc, int M>
	struct unary_expr_map<vrep_t<M>, ref_arg_holder<const_matrix<T, Mc, Nc> > >
	{
		typedef const_matrix<T, M, Nc> type;

		LMAT_ENSURE_INLINE
		static type get(const vrep_t<M>& spec,
				const ref_arg_forwarder<const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(spec.get_m(), arg.ncolumns(), arg.value());
		}
	};

	template<typename T, int Mc, int Nc, int M>
	struct unary_expr_map<vrep_t<M>, copy_arg_holder<const_matrix<T, Mc, Nc> > >
	{
		typedef const_matrix<T, M, Nc> type;

		LMAT_ENSURE_INLINE
		static type get(const vrep_t<M>& spec,
				const copy_arg_forwarder<const_matrix<T, Mc, Nc> >& arg_fwd)
		{
			typedef const_matrix<T, Mc, Nc> arg_t;
			const arg_t& arg = arg_fwd.arg;

			return type(spec.get_m(), arg.ncolumns(), arg.value());
		}
	};

	// construction functions

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<hrep_t<DynamicDim>, ref_arg_holder<Arg> >::type
	hrep(const IMatrixXpr<Arg, T>& arg, const index_t n)
	{
		return make_expr(hrep_t<DynamicDim>(n), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<hrep_t<N>, ref_arg_holder<Arg> >::type
	hrep(const IMatrixXpr<Arg, T>& arg, fixed_dim<N>)
	{
		return make_expr(hrep_t<N>(), ref_arg(arg.derived()));
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<vrep_t<DynamicDim>, ref_arg_holder<Arg> >::type
	vrep(const IMatrixXpr<Arg, T>& arg, const index_t m)
	{
		return make_expr(vrep_t<DynamicDim>(m), ref_arg(arg.derived()));
	}

	template<typename T, class Arg, int M>
	LMAT_ENSURE_INLINE
	inline typename unary_expr_map<vrep_t<M>, ref_arg_holder<Arg> >::type
	vrep(const IMatrixXpr<Arg, T>& arg, fixed_dim<M>)
	{
		return make_expr(vrep_t<M>(), ref_arg(arg.derived()));
	}


}

#endif /* REPEAT_VECS_H_ */
