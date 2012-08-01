/*
 * @file dense_matrix.h
 *
 * The matrix to represent a dense matrix.
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_MATRIX_H_
#define LIGHTMAT_DENSE_MATRIX_H_

#include "bits/dense_matrix_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols, typename Align>
	struct matrix_traits<dense_matrix<T, CTRows, CTCols, Align> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct ct_has_continuous_layout<dense_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_base_aligned<dense_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_percol_aligned<dense_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = is_same<Align, percol_aligned>::value;
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct is_linear_accessible<dense_matrix<T, CTRows, CTCols, Align> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, typename Align, class DMat>
	struct default_matrix_eval_policy<dense_matrix<T, CTRows, CTCols, Align>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	/********************************************
	 *
	 *  dense_matrix
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols, typename Align>
	class dense_matrix : public IDenseMatrix<dense_matrix<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type.");

		static_assert(CTRows >= 0 && CTCols >= 0,
				"CTRows and CTCols must be non-negative numbers.");

		static_assert(is_same<Align, base_aligned>::value,
				"Align must be base_aligned in current version.");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

	public:
		LMAT_ENSURE_INLINE dense_matrix()
		: m_internal()
		{
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n)
		: m_internal(m, n)
		{
		}

		template<class Gen>
		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n,
				const IMatrixGenerator<Gen, T>& gen)
		: m_internal(m, n)
		{
			gen.generate_to(*this);
		}

		LMAT_ENSURE_INLINE dense_matrix(const dense_matrix& s)
		: m_internal(s.m_internal)
		{
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix(const IMatrixXpr<Expr, T>& r)
		: m_internal(r.nrows(), r.ncolumns())
		{
			default_evaluate(r, *this);
		}

		LMAT_ENSURE_INLINE void swap(dense_matrix& s)
		{
			m_internal.swap(s.m_internal);
		}

	public:
		template<class Gen>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixGenerator<Gen, T>& gen)
		{
			gen.generate_to(*this);
			return *this;
		}

		LMAT_ENSURE_INLINE dense_matrix& operator = (const dense_matrix& r)
		{
			if (this != &r)
			{
				default_assign(*this, r);
			}

			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixXpr<Expr, T>& r)
		{
			default_assign(*this, r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE index_type nelems() const
		{
			return m_internal.nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_internal.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_internal.ncolumns();
		}

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_internal.nrows();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_internal.ptr_data();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_internal.ptr_data();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_col(const index_type j) const
		{
			return ptr_data() + j * lead_dim();
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_type j)
		{
			return ptr_data() + j * lead_dim();
		}

		LMAT_ENSURE_INLINE index_type offset(const index_type i, const index_type j) const
		{
			return m_internal.offset(i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(const index_type i, const index_type j)
		{
			return ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return ptr_data()[i];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return ptr_data()[i];
		}

		LMAT_ENSURE_INLINE void require_size(index_type m, index_type n)
		{
			m_internal.resize(m, n);
		}

	private:
		detail::dense_matrix_internal<T, CTRows, CTCols> m_internal;
	};


	template<typename T, int CTRows, int CTCols>
	LMAT_ENSURE_INLINE
	inline void swap(dense_matrix<T, CTRows, CTCols>& a, dense_matrix<T, CTRows, CTCols>& b)
	{
		a.swap(b);
	}


	/********************************************
	 *
	 *  derived vectors
	 *
	 ********************************************/

	template<typename T, int CTRows, typename Align>
	class dense_col : public dense_matrix<T, CTRows, 1, Align>
	{
		typedef dense_matrix<T, CTRows, 1, Align> base_mat_t;

	public:
		LMAT_ENSURE_INLINE dense_col() : base_mat_t(CTRows, 1) { }

		LMAT_ENSURE_INLINE explicit dense_col(index_t m) : base_mat_t(m, 1) { }

		template<class Gen>
		LMAT_ENSURE_INLINE
		dense_col(index_t m, const IMatrixGenerator<Gen,T>& gen) : base_mat_t(m, 1, gen) { }

		LMAT_ENSURE_INLINE dense_col(const base_mat_t& s) : base_mat_t(s) { }

		LMAT_ENSURE_INLINE dense_col(const dense_col& s) : base_mat_t(s) { }

		template<class Other>
		LMAT_ENSURE_INLINE dense_col(const IMatrixView<Other, T>& r) : base_mat_t(r) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_col(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

	public:
		template<class Gen>
		LMAT_ENSURE_INLINE
		dense_col& operator = (const IMatrixGenerator<Gen,T>& gen)
		{
			base_mat_t::operator = (gen);
			return *this;
		}

		LMAT_ENSURE_INLINE dense_col& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_col& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE
		void require_size(index_t m, index_t n)
		{
			base_mat_t::require_size(m, n);
		}

		LMAT_ENSURE_INLINE
		void require_size(index_t n)
		{
			base_mat_t::require_size(n, 1);
		}
	};


	template<typename T, int CTCols, typename Align>
	class dense_row : public dense_matrix<T, 1, CTCols, Align>
	{
		typedef dense_matrix<T, 1, CTCols, Align> base_mat_t;

	public:
		LMAT_ENSURE_INLINE dense_row() : base_mat_t(1, CTCols) { }

		LMAT_ENSURE_INLINE explicit dense_row(index_t n) : base_mat_t(1, n) { }

		template<class Gen>
		LMAT_ENSURE_INLINE
		dense_row(index_t n, const IMatrixGenerator<Gen,T>& gen) : base_mat_t(1, n, gen) { }

		LMAT_ENSURE_INLINE dense_row(const base_mat_t& s) : base_mat_t(s) { }

		LMAT_ENSURE_INLINE dense_row(const dense_row& s) : base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_row(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

	public:
		template<class Gen>
		LMAT_ENSURE_INLINE
		dense_row& operator = (const IMatrixGenerator<Gen,T>& gen)
		{
			base_mat_t::operator = (gen);
			return *this;
		}

		LMAT_ENSURE_INLINE dense_row& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_row& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE
		void require_size(index_t m, index_t n)
		{
			base_mat_t::require_size(m, n);
		}

		LMAT_ENSURE_INLINE
		void require_size(index_t n)
		{
			base_mat_t::require_size(1, n);
		}
	};


	/********************************************
	 *
	 *  evaluation
	 *
	 ********************************************/

	template<typename T, class Expr>
	LMAT_ENSURE_INLINE
	dense_matrix<T, ct_rows<Expr>::value, ct_cols<Expr>::value>
	eval(const IMatrixXpr<Expr, T>& expr)
	{
		return dense_matrix<T, ct_rows<Expr>::value, ct_cols<Expr>::value>(
				expr.derived());
	}

	template<typename T, class Expr, class Context>
	LMAT_ENSURE_INLINE
	dense_matrix<T, ct_rows<Expr>::value, ct_cols<Expr>::value>
	eval(const IMatrixXpr<Expr, T>& expr, const Context& ctx)
	{
		dense_matrix<T, ct_rows<Expr>::value, ct_cols<Expr>::value> r(expr.nrows(), expr.ncolumns());
		evaluate(expr.derived(), r, ctx);
	}

	template<typename T, class Expr>
	LMAT_ENSURE_INLINE
	T to_scalar(const IMatrixXpr<Expr, T>& expr)
	{
		dense_matrix<T,1,1> r(expr);
		return r[0];
	}


	/********************************************
	 *
	 *  Typedefs
	 *
	 ********************************************/

	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat, DynamicDim, DynamicDim)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat22, 2, 2)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat23, 2, 3)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat32, 3, 2)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat33, 3, 3)

	LMAT_MATRIX_TYPEDEFS1(dense_col, col, DynamicDim)
	LMAT_MATRIX_TYPEDEFS1(dense_col, col2, 2)
	LMAT_MATRIX_TYPEDEFS1(dense_col, col3, 3)

	LMAT_MATRIX_TYPEDEFS1(dense_row, row, DynamicDim)
	LMAT_MATRIX_TYPEDEFS1(dense_row, row2, 2)
	LMAT_MATRIX_TYPEDEFS1(dense_row, row3, 3)

}


#endif 
