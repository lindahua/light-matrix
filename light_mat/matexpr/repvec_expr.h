/**
 * @file repvec_expr.h
 *
 * @brief Vector-repeating expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REPVEC_EXPR_H_
#define LIGHTMAT_REPVEC_EXPR_H_

#include <light_mat/mateval/ewise_eval.h>

namespace lmat
{
	// forward declarations

	template<class Arg, int CN=0> class repcol_expr;
	template<class Arg, int CM=0> class reprow_expr;


	/********************************************
	 *
	 *  matrix trait classes
	 *
	 ********************************************/

	template<class Arg, int CN>
	struct matrix_traits<repcol_expr<Arg, CN> >
	: public matrix_xpr_traits_base<
	  typename meta::value_type_of<Arg>::type,
	  meta::nrows<Arg>::value,
	  CN,
	  typename meta::domain_of<Arg>::type> { };


	template<class Arg, int CM>
	struct matrix_traits<reprow_expr<Arg, CM> >
	: public matrix_xpr_traits_base<
	  typename meta::value_type_of<Arg>::type,
	  CM,
	  meta::ncols<Arg>::value,
	  typename meta::domain_of<Arg>::type> { };



	/********************************************
	 *
	 *  matrix expression classes
	 *
	 ********************************************/

	template<typename Arg, int CN>
	class repcol_expr
	: public ewise_matrix_base<repcol_expr<Arg, CN> >
	{
		static_assert( meta::is_regular_mat<Arg>::value, "Arg must be a regular matrix" );
		typedef ewise_matrix_base<repcol_expr<Arg, CN> > base_t;

	public:
		LMAT_ENSURE_INLINE
		repcol_expr(const Arg& a, dimension<CN> rdim)
		: base_t(a.nrows(), rdim.value()), m_arg(a)
		{
			LMAT_CHECK_DIMS( is_column(a) )
		}

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	template<typename Arg, int CM>
	class reprow_expr
	: public ewise_matrix_base<reprow_expr<Arg, CM> >
	{
		static_assert( meta::is_regular_mat<Arg>::value, "Arg must be a regular matrix" );
		typedef ewise_matrix_base<reprow_expr<Arg, CM> > base_t;

	public:
		LMAT_ENSURE_INLINE
		reprow_expr(const Arg& a, dimension<CM> cdim)
		: base_t(cdim.value(), a.ncolumns()), m_arg(a)
		{
			LMAT_CHECK_DIMS( is_row(a) )
		}

		LMAT_ENSURE_INLINE const Arg& arg() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	/********************************************
	 *
	 *  Expression construction functions
	 *
	 ********************************************/

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline repcol_expr<Arg, 0> repcol(const IRegularMatrix<Arg, T>& arg, index_t n)
	{
		return repcol_expr<Arg, 0>(arg.derived(), n);
	}

	template<typename T, class Arg, int CN>
	LMAT_ENSURE_INLINE
	inline repcol_expr<Arg, CN> repcol(const IRegularMatrix<Arg, T>& arg, dimension<CN> n)
	{
		return repcol_expr<Arg, CN>(arg.derived(), n);
	}

	template<typename T, class Arg>
	LMAT_ENSURE_INLINE
	inline reprow_expr<Arg, 0> reprow(const IRegularMatrix<Arg, T>& arg, index_t m)
	{
		return reprow_expr<Arg, 0>(arg.derived(), m);
	}

	template<typename T, class Arg, int CM>
	LMAT_ENSURE_INLINE
	inline reprow_expr<Arg, CM> reprow(const IRegularMatrix<Arg, T>& arg, dimension<CM> m)
	{
		return reprow_expr<Arg, CM>(arg.derived(), m);
	}


	/********************************************
	 *
	 *  Accessor classes
	 *
	 ********************************************/

	namespace internal
	{
		template<typename Arg, int CN, typename U>
		struct multicol_reader_map<repcol_expr<Arg, CN>, U>
		{
			typedef repcol_expr<Arg, CN> expr_type;
			typedef typename repeatcol_reader_map<Arg, U>::type type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return make_multicol_accessor(U(), repcol_(expr.arg()));
			}
		};

		template<typename Arg, int CM, typename U>
		struct multicol_reader_map<reprow_expr<Arg, CM>, U>
		{
			typedef reprow_expr<Arg, CM> expr_type;
			typedef typename repeatrow_reader_map<Arg, U>::type type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return make_multicol_accessor(U(), reprow_(expr.arg()));
			}
		};
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	template<typename Arg, int CN>
	struct supports_linear_access<repcol_expr<Arg, CN> >
	{
		static const bool value = false;
	};

	template<typename Arg, int CM>
	struct supports_linear_access<reprow_expr<Arg, CM> >
	{
		static const bool value = false;
	};

	template<typename Arg, int CN, typename Kind>
	struct supports_simd<repcol_expr<Arg, CN>, Kind>
	: public supports_simd<Arg, Kind> { };

	template<typename Arg, int CM, typename Kind>
	struct supports_simd<reprow_expr<Arg, CM>, Kind>
	: public supports_simd<typename matrix_traits<Arg>::value_type, Kind> { };


	template<typename Arg, int CN, class DMat>
	inline void evaluate(const repcol_expr<Arg, CN>& sexpr,
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


	template<typename Arg, int CM, class DMat>
	inline void evaluate(const reprow_expr<Arg, CM>& sexpr,
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
