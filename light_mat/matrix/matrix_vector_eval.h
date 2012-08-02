/**
 * @file matrix_vector_eval.h
 *
 * Generic vector-based matrix evaluation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_VECTOR_EVAL_H_
#define LIGHTMAT_MATRIX_VECTOR_EVAL_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/matrix/matrix_classes.h>

namespace lmat
{
	// Policy types

	struct by_scalars { };
	struct by_simd { };

	struct as_linear_vec { };
	struct per_column { };

	// Means can be either by_scalars or by_simd
	// Org can be either as_linear_vec or per_column

	template<typename Org, typename Means> struct vector_eval_policy { };
	template<class Expr, typename Org, typename Means> struct vector_eval;

	template<class Expr>
	struct vector_eval_default_policy
	{
		// Note: SIMD has not been implemented

		typedef by_scalars means;

		static const int linear_cost = vector_eval<Expr, as_linear_vec, means>::cost;
		static const int percol_cost = vector_eval<Expr, per_column, means>::cost;

		static const bool choose_linear = (linear_cost <= percol_cost);

		typedef typename if_c<choose_linear, as_linear_vec, per_column>::type org;

		typedef vector_eval_policy<org, means> type;
	};


	/********************************************
	 *
	 *  Vector evaluators
	 *
	 ********************************************/

	// interfaces


	template<class Derived, typename T>
	class ILinearVectorEvaluator
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_value(const index_t i) const
		{
			return derived().get_value(i);
		}
	};


	template<class Evaluator> struct percol_eval_state;

	struct nil_eval_state { };

	template<class Derived, typename T>
	class IPerColVectorEvaluator
	{
	public:
		LMAT_CRTP_REF

		LMAT_ENSURE_INLINE
		T get_value(const typename percol_eval_state<Derived>::type& col_stat,
				const index_t i) const
		{
			return derived().get_value(col_stat, i);
		}

		LMAT_ENSURE_INLINE
		typename percol_eval_state<Derived>::type col_state(const index_t j) const
		{
			return derived().col_state(j);
		}
	};


	template<class Eva, typename T>
	struct is_linear_vector_evaluator
	{
		static const bool value = is_base_of<
				ILinearVectorEvaluator<Eva, T>,
				Eva>::value;
	};

	template<class Eva, typename T>
	struct is_percol_vector_evaluator
	{
		static const bool value = is_base_of<
				IPerColVectorEvaluator<Eva, T>,
				Eva>::value;
	};


	// forward declarations

	template<typename T> class continuous_linear_evaluator;
	template<typename T> class dense_percol_evaluator;
	template<typename T> class const_linear_evaluator;
	template<typename T> class const_percol_evaluator;
	template<typename T> class cached_linear_evaluator;
	template<typename T> class cached_percol_evaluator;

	// specific evaluators

	template<typename T>
	class continuous_linear_evaluator
	: public ILinearVectorEvaluator<continuous_linear_evaluator<T>, T>
	{
	public:

		template<class Mat>
		LMAT_ENSURE_INLINE
		continuous_linear_evaluator(const IDenseMatrix<Mat, T>& X)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(ct_has_continuous_layout<Mat>::value,
					"Mat must always have continuous layout");
#endif
			m_data = X.ptr_data();
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}
	private:
		const T *m_data;
	};


	template<typename T>
	class densecol_state
	{
	private:
		const T* m_ptrcol;

	public:
		LMAT_ENSURE_INLINE
		densecol_state(const T* p) : m_ptrcol(p) { }

		LMAT_ENSURE_INLINE
		const T *ptr() const { return m_ptrcol; }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t i) const
		{
			return m_ptrcol[i];
		}
	};

	template<typename T>
	struct percol_eval_state<dense_percol_evaluator<T> >
	{
		typedef densecol_state<T> type;
	};

	template<typename T>
	class dense_percol_evaluator
	: public IPerColVectorEvaluator<dense_percol_evaluator<T>, T>
	{
	public:
		template<class Mat>
		LMAT_ENSURE_INLINE
		dense_percol_evaluator(const IDenseMatrix<Mat, T>& X)
		: m_ldim(X.lead_dim())
		, m_data(X.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const densecol_state<T>& s, const index_t i) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE densecol_state<T> col_state(const index_t j) const
		{
			return m_data + m_ldim * j;
		}
	private:
		const index_t m_ldim;
		const T *m_data;
	};


	template<typename T>
	class const_linear_evaluator
	: public ILinearVectorEvaluator<const_linear_evaluator<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_linear_evaluator(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_val;
		}

	private:
		const T m_val;
	};

	template<typename T>
	struct percol_eval_state<const_percol_evaluator<T> >
	{
		typedef nil_eval_state type;
	};

	template<typename T>
	class const_percol_evaluator
	: public IPerColVectorEvaluator<const_percol_evaluator<T>, T>
	{
	public:
		template<int CTRows, int CTCols>
		LMAT_ENSURE_INLINE
		const_percol_evaluator(const const_matrix<T, CTRows, CTCols>& X)
		: m_val(X.value())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const nil_eval_state&, const index_t i) const
		{
			return m_val;
		}

		LMAT_ENSURE_INLINE nil_eval_state col_state(const index_t ) const
		{
			return nil_eval_state();
		}

	private:
		const T m_val;
	};



	template<typename T>
	class cached_linear_evaluator
	: public ILinearVectorEvaluator<cached_linear_evaluator<T>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_linear_evaluator(const IMatrixXpr<Expr, T>& X)
		: m_cache(X), m_data(m_cache.ptr_data())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const index_t i) const
		{
			return m_data[i];
		}
	private:
		dense_matrix<T> m_cache;
		const T *m_data;
	};


	template<typename T>
	struct percol_eval_state<cached_percol_evaluator<T> >
	{
		typedef densecol_state<T> type;
	};

	template<typename T>
	class cached_percol_evaluator
	: public IPerColVectorEvaluator<cached_percol_evaluator<T>, T>
	{
	public:
		template<class Expr>
		LMAT_ENSURE_INLINE
		cached_percol_evaluator(const IMatrixXpr<Expr, T>& X)
		: m_cache(X), m_ldim(m_cache.lead_dim())
		{
		}

		LMAT_ENSURE_INLINE T get_value(const densecol_state<T>& s, const index_t i) const
		{
			return s[i];
		}

		LMAT_ENSURE_INLINE densecol_state<T> col_state(const index_t j) const
		{
			return m_cache.ptr_data() + m_ldim * j;
		}

	private:
		dense_matrix<T> m_cache;
		const index_t m_ldim;
	};


	/********************************************
	 *
	 *  evaluation context
	 *
	 ********************************************/

	template<typename T, int CTSize, typename Means> struct linear_eval_impl;
	template<typename T, int CTRows, int CTCols, typename Means> struct percol_eval_impl;

	template<class Expr, class Dst, typename Org, typename Means> struct vector_eval_impl_map;

	template<class Expr, class Dst, typename Means>
	struct vector_eval_impl_map<Expr, Dst, as_linear_vec, Means>
	{
		typedef linear_eval_impl<
			typename matrix_traits<Expr>::value_type,
			binary_ct_size<Expr, Dst>::value,
			Means> type;
	};

	template<class Expr, class Dst, typename Means>
	struct vector_eval_impl_map<Expr, Dst, per_column, Means>
	{
		typedef percol_eval_impl<
			typename matrix_traits<Expr>::value_type,
			binary_ct_rows<Expr, Dst>::value,
			binary_ct_cols<Expr, Dst>::value,
			Means> type;
	};



	template<typename T, int CTSize>
	struct linear_eval_impl<T, CTSize, by_scalars>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTSize > 0, "CTSize must be positive.");
#endif

		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			T* pd = dst.ptr_data();
			for (index_t i = 0; i < CTSize; ++i)
			{
				pd[i] = evaluator.get_value(i);
			}
		}
	};

	template<typename T>
	struct linear_eval_impl<T, DynamicDim, by_scalars>
	{
		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				const ILinearVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			T* pd = dst.ptr_data();
			const index_t len = dst.nelems();

			for (index_t i = 0; i < len; ++i)
			{
				pd[i] = evaluator.get_value(i);
			}
		}
	};


	template<typename T, int CTRows, int CTCols>
	struct percol_eval_impl<T, CTRows, CTCols, by_scalars>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(CTRows > 0, "CTSize must be positive.");
#endif

		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				typename percol_eval_state<Eva>::type s = evaluator.col_state(j);

				for (index_t i = 0; i < CTRows; ++i)
				{
					pd[i] = evaluator.get_value(s, i);
				}
			}
		}
	};


	template<typename T, int CTCols>
	struct percol_eval_impl<T, DynamicDim, CTCols, by_scalars>
	{
		template<class Eva, class Mat>
		LMAT_ENSURE_INLINE
		static void evaluate(
				IPerColVectorEvaluator<Eva, T>& evaluator,
				IDenseMatrix<Mat, T>& dst)
		{
			const index_t nrows = dst.nrows();
			const index_t ncols = dst.ncolumns();
			const index_t ldim = dst.lead_dim();
			T *pd = dst.ptr_data();

			for (index_t j = 0; j < ncols; ++j, pd += ldim)
			{
				typename percol_eval_state<Eva>::type s = evaluator.col_state(j);

				for (index_t i = 0; i < nrows; ++i)
				{
					pd[i] = evaluator.get_value(s, i);
				}
			}
		}
	};


	/********************************************
	 *
	 *  Evaluator Map and Cost Model
	 *
	 ********************************************/

	const int VEC_EVAL_CACHE_COST = 1000;
	const int SHORTVEC_LENGTH_THRESHOLD = 4;
	const int SHORTVEC_PERCOL_COST = 200;

	template<class Expr, typename Org, typename Means> struct generic_vector_eval;

	template<class Expr>
	struct generic_vector_eval<Expr, as_linear_vec, by_scalars>
	{
		static const bool can_direct = is_dense_mat<Expr>::value && ct_has_continuous_layout<Expr>::value;

		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_c<can_direct,
					continuous_linear_evaluator<T>,
					cached_linear_evaluator<T> >::type evaluator_type;

		static const int cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
	};

	template<class Expr>
	struct generic_vector_eval<Expr, per_column, by_scalars>
	{
		static const bool can_direct = is_dense_mat<Expr>::value;
		static const bool has_short_col = ct_rows<Expr>::value < SHORTVEC_LENGTH_THRESHOLD;

		typedef typename matrix_traits<Expr>::value_type T;

		typedef typename
				if_c<can_direct,
					dense_percol_evaluator<T>,
					cached_percol_evaluator<T> >::type evaluator_type;

		static const int normal_cost = can_direct ? 0 : VEC_EVAL_CACHE_COST;
		static const int shortv_cost = SHORTVEC_PERCOL_COST + normal_cost;

		static const int cost = has_short_col ? shortv_cost : normal_cost;
	};

	template<class Expr, typename Means>
	struct vector_eval<Expr, as_linear_vec, Means>
	{
		typedef typename generic_vector_eval<Expr, as_linear_vec, Means>::evaluator_type evaluator_type;
		static const int cost = generic_vector_eval<Expr, as_linear_vec, Means>::cost;
	};

	template<class Expr, typename Means>
	struct vector_eval<Expr, per_column, Means>
	{
		typedef typename generic_vector_eval<Expr, per_column, Means>::evaluator_type evaluator_type;
		static const int normal_cost = generic_vector_eval<Expr, per_column, Means>::normal_cost;
		static const int shortv_cost = generic_vector_eval<Expr, per_column, Means>::shortv_cost;
		static const int cost = generic_vector_eval<Expr, per_column, Means>::cost;
	};


	template<typename T, int CTRows, int CTCols>
	struct vector_eval<const_matrix<T, CTRows, CTCols>, as_linear_vec, by_scalars>
	{
		typedef const_linear_evaluator<T> evaluator_type;

		static const int cost = 0;
	};

	template<typename T, int CTRows, int CTCols>
	struct vector_eval<const_matrix<T, CTRows, CTCols>, per_column, by_scalars>
	{
		typedef const_percol_evaluator<T> evaluator_type;

		static const int normal_cost = 0;
		static const int shortv_cost = 0;
		static const int cost = 0;
	};


	/********************************************
	 *
	 *  Evaluation functions
	 *
	 ********************************************/

	template<typename T, class Expr, class Dst, typename Org, typename Means>
	LMAT_ENSURE_INLINE
	void evaluate(const IMatrixXpr<Expr, T>& src, IDenseMatrix<Dst, T>& dst, vector_eval_policy<Org, Means>)
	{
		typedef typename vector_eval_impl_map<Expr, Dst, Org, Means>::type impl_t;
		typename vector_eval<Expr, Org, Means>::evaluator_type evaluator(src.derived());
		impl_t::evaluate(evaluator, dst.derived());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void linear_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), vector_eval_policy<as_linear_vec, by_scalars>());
	}

	template<typename T, class Expr, class Dst>
	LMAT_ENSURE_INLINE
	void percol_by_scalars_evaluate(const IMatrixXpr<Expr, T>& expr, IDenseMatrix<Dst, T>& dst)
	{
		evaluate(expr.derived(), dst.derived(), vector_eval_policy<per_column, by_scalars>());
	}

}

#endif
