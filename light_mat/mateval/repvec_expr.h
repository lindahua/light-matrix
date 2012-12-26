/**
 * @file repvec_expr.h
 *
 * @brief Vector-repeating expressions
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_REPVEC_EXPR_H_
#define LIGHTMAT_REPVEC_EXPR_H_

#include "internal/map_expr_internal.h"

namespace lmat
{
	// forward declarations

	template<class Arg, int N=0> class repcol_expr;
	template<class Arg, int M=0> class reprow_expr;


	/********************************************
	 *
	 *  matrix trait classes
	 *
	 ********************************************/

	template<class Arg, int N>
	struct matrix_traits<repcol_expr<Arg, N> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = meta::nrows<Arg>::value;
		static const int ct_num_cols = N;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};


	template<class Arg, int M>
	struct matrix_traits<reprow_expr<Arg, M> >
	{
		static const int num_dimensions = 2;
		static const int ct_num_rows = M;
		static const int ct_num_cols = meta::ncols<Arg>::value;

		static const bool is_readonly = true;

		typedef matrix_shape<ct_num_rows, ct_num_cols> shape_type;

		typedef typename matrix_traits<Arg>::value_type value_type;
		typedef typename matrix_traits<Arg>::domain domain;
	};



	/********************************************
	 *
	 *  matrix expression classes
	 *
	 ********************************************/

	template<typename Arg, int N>
	class repcol_expr
	: public IEWiseMatrix<repcol_expr<Arg, N>, typename matrix_traits<Arg>::value_type>
	{
		static_assert( meta::is_regular_mat<Arg>::value, "Arg must be a regular matrix" );

		typedef matrix_shape<meta::nrows<Arg>::value, N> shape_type;

	public:
		LMAT_ENSURE_INLINE
		repcol_expr(const Arg& a, dimension<N> rdim)
		: m_shape(a.nrows(), rdim.value()), m_arg(a)
		{
			LMAT_CHECK_DIMS( is_column(a) )
		}

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
		const Arg& m_arg;
	};


	template<typename Arg, int M>
	class reprow_expr
	: public IEWiseMatrix<reprow_expr<Arg, M>, typename matrix_traits<Arg>::value_type>
	{
		static_assert( meta::is_regular_mat<Arg>::value, "Arg must be a regular matrix" );

		typedef matrix_shape<M, meta::ncols<Arg>::value> shape_type;

	public:
		LMAT_ENSURE_INLINE
		reprow_expr(const Arg& a, dimension<M> cdim)
		: m_shape(cdim.value(), a.ncolumns()), m_arg(a)
		{
			LMAT_CHECK_DIMS( is_row(a) )
		}

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
		const Arg& m_arg;
	};


	/********************************************
	 *
	 *  Expression construction functions
	 *
	 ********************************************/

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	repcol_expr<Arg, 0> repcol(const IRegularMatrix<Arg, T>& arg, index_t n)
	{
		return repcol_expr<Arg, 0>(arg.derived(), n);
	}

	template<typename T, class Arg, int N>
	LMAT_ENSURE_INLINE
	repcol_expr<Arg, N> repcol(const IRegularMatrix<Arg, T>& arg, dimension<N> n)
	{
		return repcol_expr<Arg, N>(arg.derived(), n);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	reprow_expr<Arg, 0> reprow(const IRegularMatrix<Arg, T>& arg, index_t m)
	{
		return reprow_expr<Arg, 0>(arg.derived(), m);
	}

	template<typename T, class Arg, int M>
	LMAT_ENSURE_INLINE
	reprow_expr<Arg, M> reprow(const IRegularMatrix<Arg, T>& arg, dimension<M> m)
	{
		return reprow_expr<Arg, M>(arg.derived(), m);
	}


	/********************************************
	 *
	 *  Accessor classes
	 *
	 ********************************************/

	namespace internal
	{
		template<typename Arg, int N, typename U>
		struct multicol_reader_map<repcol_expr<Arg, N>, U>
		{
			typedef repcol_expr<Arg, N> expr_type;
			typedef typename repeatcol_reader_map<Arg, U>::type type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return make_multicol_accessor(U(), in_(expr.arg(), atags::repcol()));
			}
		};

		template<typename Arg, int M, typename U>
		struct multicol_reader_map<reprow_expr<Arg, M>, U>
		{
			typedef reprow_expr<Arg, M> expr_type;
			typedef typename repeatrow_reader_map<Arg, U>::type type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return make_multicol_accessor(U(), in_(expr.arg(), atags::reprow()));
			}
		};
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename Arg, int N>
		struct prefers_linear<repcol_expr<Arg, N> >
		{
			static const bool value = false;
		};

		template<typename Arg, int M>
		struct prefers_linear<reprow_expr<Arg, M> >
		{
			static const bool value = false;
		};

		template<typename Arg, int N,
			typename T, typename Kind, bool IsLinear>
		struct prefers_simd<repcol_expr<Arg, N>, T, Kind, IsLinear>
		{
			static const bool value = arg_prefers_simd<Arg, Kind, IsLinear>::value;
		};

		template<typename Arg, int M,
			typename T, typename Kind, bool IsLinear>
		struct prefers_simd<reprow_expr<Arg, M>, T, Kind, IsLinear>
		{
			static const bool value = arg_prefers_simd<Arg, Kind, IsLinear>::value;
		};
	}

	template<typename Arg, int N, class DMat>
	inline void evaluate(const repcol_expr<Arg, N>& sexpr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		const Arg& a = sexpr.arg();
		const index_t m = sexpr.nrows();
		const index_t n = sexpr.ncolumns();

		if (n == 1)
		{
			auto dcol = dmat.column(0);
			copy(a, dcol);
		}
		else if (m == 1)
		{
			fill(dmat, *(a.ptr_data()));
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				auto dcol = dmat.column(j);
				copy(a, dcol);
			}
		}
	}


	template<typename Arg, int M, class DMat>
	inline void evaluate(const reprow_expr<Arg, M>& sexpr,
			IRegularMatrix<DMat, typename matrix_traits<Arg>::value_type>& dmat)
	{
		typedef typename matrix_traits<Arg>::value_type T;

		const Arg& a = sexpr.arg();
		const index_t m = sexpr.nrows();
		const index_t n = sexpr.ncolumns();

		if (m == 1)
		{
			auto drow = dmat.row(0);
			copy(a, drow);
		}
		else if (n == 1)
		{
			fill(dmat, *(a.ptr_data()));
		}
		else
		{
			for (index_t j = 0; j < n; ++j)
			{
				auto dcol = dmat.column(j);
				fill(dcol, a(0, j));
			}
		}
	}


}

#endif /* REPVEC_EXPR_H_ */
