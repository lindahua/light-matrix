/**
 * @file matrix_access_base.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_ACCESS_BASE_H_
#define LIGHTMAT_MATRIX_ACCESS_BASE_H_

#include "bits/macc_eval_impl.h"

namespace lmat
{

	/********************************************
	 *
	 *  Accessor Interfaces
	 *
	 ********************************************/

	template<class SExpr, class AccCate, class KerCate>
	struct macc_accessor_map;

	template<class Derived, typename T>
	class ILinearMatrixScalarAccessor
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return derived().get_scalar(i);
		}
	};


	template<class Accessor> struct percol_macc_state_map;

	template<class Derived, typename T>
	class IPerColMatrixScalarAccessor
	{
	public:
		typedef typename percol_macc_state_map<Derived>::type col_state_t;

		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, const col_state_t& s) const
		{
			return derived().get_scalar(i, s);
		}

		LMAT_ENSURE_INLINE
		col_state_t col_state(const index_t j) const
		{
			return derived().col_state(j);
		}
	};


	/********************************************
	 *
	 *  Accessors for scalar expression
	 *
	 ********************************************/

	template<typename T> class const_scalar_accessor;

	// scalar accessing

	template<class T>
	struct percol_macc_state_map<const_scalar_accessor<T> >
	{
		typedef nil_t type;
	};

	template<typename T>
	class const_scalar_accessor
		: public ILinearMatrixScalarAccessor<const_scalar_accessor<T>, T>
		, public IPerColMatrixScalarAccessor<const_scalar_accessor<T>, T>
	{
	public:
		LMAT_ENSURE_INLINE
		const_scalar_accessor(const scalar_expr<T>& xpr)
		: m_val(xpr.value)
		{ }

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i, nil_t) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE
		nil_t col_state(const index_t j) const
		{
			return nil_t();
		}

	private:
		const T m_val;
	};


	template<typename T, class AccCate>
	struct macc_accessor_map<scalar_expr<T>, AccCate, scalar_kernel_t>
	{
		typedef const_scalar_accessor<T> type;
	};


	/********************************************
	 *
	 *  macc schemes
	 *
	 ********************************************/

	template<class Xpr, typename AccCate, typename KerCate>
	struct macc_cost;

	const int MACC_CACHE_COST = 1200;
	const int MACC_SHORT_PERCOL_COST = 100;
	const int MACC_SHORTCOL_UBOUND = 4;

	// matrix access categories

	struct any_macc { };
	struct linear_macc { };
	struct percol_macc { };

	// matrix access setting

	template<typename AccCate, typename KerCate, int M, int N>
	struct macc_scheme;

	template<int M, int N>
	struct macc_scheme<any_macc, any_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_linear;
		const bool _use_simd;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool lin, bool simd)
		: shape(m, n), _use_linear(lin), _use_simd(simd) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return _use_linear; }
		LMAT_ENSURE_INLINE bool use_simd() const { return _use_simd; }
	};

	template<int M, int N>
	struct macc_scheme<linear_macc, any_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_simd;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool simd)
		: shape(m, n), _use_simd(simd) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return true; }
		LMAT_ENSURE_INLINE bool use_simd() const { return _use_simd; }
	};

	template<int M, int N>
	struct macc_scheme<percol_macc, any_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_simd;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool simd)
		: shape(m, n), _use_simd(simd) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return false; }
		LMAT_ENSURE_INLINE bool use_simd() const { return _use_simd; }
	};

	template<int M, int N>
	struct macc_scheme<any_macc, scalar_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_linear;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool lin)
		: shape(m, n), _use_linear(lin) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return _use_linear; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			if (_use_linear)
			{
				internal::macc_eval_linear_scalar::evaluate(shape.nelems(), sexpr, dmat);
			}
			else
			{
				internal::macc_eval_percol_scalar::evaluate(shape.nrows(), shape.ncolumns(), sexpr, dmat);
			}
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			const index_t m = dmat.nrows();

			int macc_linear_cost = macc_cost<SExpr, linear_macc, scalar_kernel_t>::value;
			int macc_percol_cost = macc_cost<SExpr, percol_macc, scalar_kernel_t>::value;

			if (m <= MACC_SHORT_PERCOL_COST)
				macc_percol_cost += MACC_SHORT_PERCOL_COST;

			return macc_scheme(m, dmat.ncolumns(),
					macc_linear_cost <= macc_percol_cost);
		}
	};

	template<int M, int N>
	struct macc_scheme<linear_macc, scalar_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return true; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			internal::macc_eval_linear_scalar::evaluate(shape.nelems(), sexpr, dmat);
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			return macc_scheme(dmat.nrows(), dmat.ncolumns());
		}
	};


	template<int M, int N>
	struct macc_scheme<percol_macc, scalar_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return false; }
		LMAT_ENSURE_INLINE bool use_simd() const { return false; }

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		void evaluate(const SExpr& sexpr, DMat& dmat)
		{
			internal::macc_eval_percol_scalar::evaluate(shape.nrows(), shape.ncolumns(), sexpr, dmat);
		}

		template<class SExpr, class DMat>
		LMAT_ENSURE_INLINE
		static macc_scheme get_default(const SExpr& expr, const DMat& dmat)
		{
			return macc_scheme(dmat.nrows(), dmat.ncolumns());
		}
	};


	template<int M, int N>
	struct macc_scheme<any_macc, simd_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;
		const bool _use_linear;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n, bool lin)
		: shape(m, n), _use_linear(lin) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return _use_linear; }
		LMAT_ENSURE_INLINE bool use_simd() const { return true; }
	};

	template<int M, int N>
	struct macc_scheme<linear_macc, simd_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return true; }
		LMAT_ENSURE_INLINE bool use_simd() const { return true; }
	};

	template<int M, int N>
	struct macc_scheme<percol_macc, simd_kernel_t, M, N>
	{
		const matrix_shape<M, N> shape;

		LMAT_ENSURE_INLINE
		macc_scheme(index_t m, index_t n)
		: shape(m, n) { }

		LMAT_ENSURE_INLINE bool use_linear() const { return false; }
		LMAT_ENSURE_INLINE bool use_simd() const { return true; }
	};


	template<class SExpr, class DMat>
	struct default_macc_scheme
	{
		typedef typename
				if_<ct_supports_linear_index<DMat>,
				any_macc,
				percol_macc>::type access_category;

		typedef scalar_kernel_t kernel_category;

		static const int M = common_ctrows<SExpr, DMat>::value;
		static const int N = common_ctcols<SExpr, DMat>::value;

		typedef macc_scheme<access_category, kernel_category, M, N> type;

		static type get(const SExpr& sexpr, const DMat& dmat)
		{
			return type::get_default(sexpr, dmat);
		}
	};

}

#endif /* MATRIX_ACCESS_BASE_H_ */





