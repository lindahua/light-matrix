/**
 * @file matrix_visitors.h
 *
 * The adapting classes to visit matrix elements
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VISITORS_H_
#define LIGHTMAT_MATRIX_VISITORS_H_

#include <light_mat/matrix/matrix_classes.h>

namespace lmat
{

	/********************************************
	 *
	 *  Auxiliary types
	 *
	 ********************************************/


	// kernel categories

	struct scalar_kernel_t { };
	struct simd_kernel_t { };

	// visitor categories

	struct linear_vis;
	struct percol_vis;

	// maps

	template<class Visitor>
	struct matrix_visitor_state { };

	template<typename VisCate, typename KerCate>
	struct matrix_visit_policy
	{
		typedef VisCate visitor_category;
		typedef KerCate kernel_category;
	};

	template<typename VisCate, typename KerCate, int ME, int NE>
	struct matrix_visit_setting
	{
		typedef matrix_visit_policy<VisCate, KerCate> policy_type;
		typedef VisCate visitor_category;
		typedef KerCate kernel_category;

		 static const int ct_rows = ME;
		 static const int ct_cols = NE;
	};


	/********************************************
	 *
	 *  Interfaces
	 *
	 ********************************************/

	template<class Derived, typename T>
	class ILinearMatrixScalarVisitor
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return derived().get_scalar(i);
		}
	};


	template<class Derived, typename T>
	class IPerColMatrixScalarVisitor
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i,
				const typename matrix_visitor_state<Derived>::type& s) const
		{
			return derived().get_scalar(i, s);
		}

		LMAT_ENSURE_INLINE
		typename matrix_visitor_state<Derived>::type
		col_state(const index_t j) const
		{
			return derived().col_state(j);
		}
	};


	/********************************************
	 *
	 *  Useful visitors
	 *
	 ********************************************/

	// forward declarations

	template<typename T> class continuous_linear_mvisitor;
	template<typename T> class dense_percol_mvisitor;
	template<typename T> class const_linear_mvisitor;
	template<typename T> class const_percol_mvisitor;
	template<typename T, int M, int N> class cached_linear_mvisitor;
	template<typename T, int M, int N> class cached_percol_mvisitor;

	// class definitions

	// continuous linear

	template<typename T>
	class continuous_linear_mvisitor
	: public ILinearMatrixScalarVisitor<continuous_linear_mvisitor<T>, T>
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		continuous_linear_mvisitor(const IDenseMatrix<Mat, T>& X)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(ct_has_continuous_layout<Mat>::value,
					"Mat must always have continuous layout");
#endif
			m_data = X.ptr_data();
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_data[i];
		}

	private:
		const T *m_data;
	};


	// dense percol

	template<typename T>
	struct matrix_visitor_state<dense_percol_mvisitor<T> >
	{
		typedef const T* type;
	};

	template<typename T>
	class dense_percol_mvisitor
	: public IPerColMatrixScalarVisitor<dense_percol_mvisitor<T>, T>
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		dense_percol_mvisitor(const IDenseMatrix<Mat, T>& X)
		: m_ldim(X.lead_dim())
		, m_data(X.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_scalar(const index_t i, const T* s) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE const T* col_state(const index_t j) const
		{
			return m_data + m_ldim * j;
		}
	private:
		const index_t m_ldim;
		const T *m_data;
	};


	// const linear

	template<typename T>
	class const_linear_mvisitor
	: public ILinearMatrixScalarVisitor<const_linear_mvisitor<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_linear_mvisitor(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_val;
		}
	private:
		const T m_val;
	};


	// const percol

	template<typename T>
	struct matrix_visitor_state<const_percol_mvisitor<T> >
	{
		typedef nil_type type;
	};

	template<typename T>
	class const_percol_mvisitor
	: public IPerColMatrixScalarVisitor<const_percol_mvisitor<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_percol_mvisitor(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE T get_scalar(const index_t, nil_type) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE nil_type col_state(const index_t) const
		{
			return nil_type();
		}

	private:
		const T m_val;
	};


	// cached linear

	template<typename T, int M, int N>
	class cached_linear_mvisitor
	: public ILinearMatrixScalarVisitor<cached_linear_mvisitor<T, M, N>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_linear_mvisitor(const IMatrixXpr<Expr, T>& X)
		: m_cache(X), m_data(m_cache.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_data[i];
		}
	private:
		dense_matrix<T, M, N> m_cache;
		const T* m_data;
	};


	// cached percol

	template<typename T, int M, int N>
	struct matrix_visitor_state<cached_percol_mvisitor<T, M, N> >
	{
		typedef const T* type;
	};

	template<typename T, int M, int N>
	class cached_percol_mvisitor
	: public IPerColMatrixScalarVisitor<cached_percol_mvisitor<T, M, N>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_percol_mvisitor(const IMatrixXpr<Expr, T>& X)
		: m_cache(X)
		{
		}

		LMAT_ENSURE_INLINE T get_scalar(const index_t i, const T* s) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE const T* col_state(const index_t j) const
		{
			return m_cache.ptr_col(j);
		}

	private:
		dense_matrix<T, M, N> m_cache;
	};


	/********************************************
	 *
	 *  Cost model and Visitor maps
	 *
	 ********************************************/

	// generic vector evaluation map

	const int MVISIT_CACHE_COST = 1000;
	const int SHORTVEC_LENGTH_THRESHOLD = 4;
	const int SHORTVEC_PERCOL_COST = 200;

	namespace detail
	{
		template<class Expr, typename VisCate, typename KerCate, int ME, int NE>
		struct _default_vismap;

		// linear

		template<class Expr, int ME, int NE>
		struct _default_vismap<Expr, linear_vis, scalar_kernel_t, ME, NE>
		{
			typedef typename matrix_traits<Expr>::value_type T;
			typedef cached_linear_mvisitor<T, ct_rows<Expr>::value, ct_cols<Expr>::value> cached_vis_type;

			typedef typename
					if_<is_linear_memory_accessible<Expr>,
						continuous_linear_mvisitor<T>,
						cached_vis_type>::type visitor_type;

			static const int cost = is_linear_memory_accessible<Expr>::value ? 0 : MVISIT_CACHE_COST;
		};


		// per-column

		template<class Expr, int ME, int NE>
		struct _default_vismap<Expr, percol_vis, scalar_kernel_t, ME, NE>
		{
			typedef typename matrix_traits<Expr>::value_type T;
			typedef cached_percol_mvisitor<T, ct_rows<Expr>::value, ct_cols<Expr>::value> cached_vis_type;

			typedef typename
					if_<is_dense_mat<Expr>,
						dense_percol_mvisitor<T>,
						cached_vis_type>::type visitor_type;

			static const int normal_cost = is_dense_mat<Expr>::value ? 0 : MVISIT_CACHE_COST;
			static const int shortv_cost = SHORTVEC_PERCOL_COST + normal_cost;

			static const int cost = ME < SHORTVEC_LENGTH_THRESHOLD ? shortv_cost : normal_cost;
		};
	}


	// generic maps

	template<class Expr, typename VisSetting>
	struct matrix_vismap
	{
		typedef detail::_default_vismap<Expr,
				typename VisSetting::visitor_category,
				typename VisSetting::kernel_category,
				VisSetting::ct_rows,
				VisSetting::ct_cols> _map_t;

		typedef typename _map_t::visitor_type type;

		static const int cost = _map_t::cost;
	};


	// specialized map for const_matrix

	template<typename T, int CTRows, int CTCols, int ME, int NE>
	struct matrix_vismap<
		const_matrix<T, CTRows, CTCols>,
		matrix_visit_setting<linear_vis, scalar_kernel_t, ME, NE> >
	{
		typedef const_linear_mvisitor<T> type;
		static const int cost = 0;
	};

	template<typename T, int CTRows, int CTCols, int ME, int NE>
	struct matrix_vismap<
		const_matrix<T, CTRows, CTCols>,
		matrix_visit_setting<percol_vis, scalar_kernel_t, ME, NE> >
	{
		typedef const_percol_mvisitor<T> type;
		static const int cost = 0;
	};


	// default policy

	template<class Expr, class Dst>
	struct default_matrix_visit_policy
	{
		// Note: SIMD has not been implemented

		static const int ME = binary_ct_rows<Expr, Dst>::value;
		static const int NE = binary_ct_cols<Expr, Dst>::value;

		typedef scalar_kernel_t kernel_t;

		static const int linear_cost = matrix_vismap<Expr,
				matrix_visit_setting<linear_vis, kernel_t, ME, NE> >::cost;

		static const int percol_cost = matrix_vismap<Expr,
				matrix_visit_setting<percol_vis, kernel_t, ME, NE> >::cost;

		static const bool choose_linear = ct_has_continuous_layout<Dst>::value && (linear_cost <= percol_cost);

		typedef typename if_c<choose_linear, linear_vis, percol_vis>::type vis_cate;

		typedef matrix_visit_policy<vis_cate, kernel_t> type;
	};


}

#endif /* MATRIX_VISITORS_H_ */
