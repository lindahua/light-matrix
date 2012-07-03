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
	 *  dense_matrix
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
	};

	template<typename T, int CTRows, int CTCols, typename Align>
	struct has_continuous_layout<dense_matrix<T, CTRows, CTCols, Align> >
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


	template<typename T, int CTRows, int CTCols, typename Align>
	class dense_matrix : public IDenseMatrix<dense_matrix<T, CTRows, CTCols, Align>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type.");

		static_assert(CTRows >= 0 && CTCols >= 0,
				"CTRows and CTCols must be non-negative numbers.");

		static_assert(is_same<Align, base_aligned>::value || is_same<Align, percol_aligned>::value,
				"Align must be either base_aligned or percol_aligned.");
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
			gen.generate_to(m, n, m, m_internal.ptr_data());
		}

		LMAT_ENSURE_INLINE dense_matrix(const dense_matrix& s)
		: m_internal(s.m_internal)
		{
		}

		template<class Other>
		LMAT_ENSURE_INLINE dense_matrix(const IDenseMatrix<Other, T>& r)
		: m_internal(r.nrows(), r.ncolumns())
		{
			copy(r.derived(), *this);
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix(const IMatrixXpr<Expr, T>& r)
		: m_internal(r.nrows(), r.ncolumns())
		{
			evaluate_to(r.derived(), *this);
		}

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_matrix(const IMatrixXpr<Expr, S>& r)
		: m_internal(r.nrows(), r.ncolumns())
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(is_implicitly_convertible<S, T>::value,
					"S is NOT implicitly-convertible to T.");
#endif
			implicitly_cast_to(r.derived(), *this);
		}

		LMAT_ENSURE_INLINE void swap(dense_matrix& s)
		{
			m_internal.swap(s.m_internal);
		}

	public:
		template<class Gen>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixGenerator<Gen, T>& gen)
		{
			gen.generate_to(nrows(), ncolumns(), lead_dim(), ptr_data());
			return *this;
		}

		LMAT_ENSURE_INLINE dense_matrix& operator = (const dense_matrix& r)
		{
			if (this != &r)
			{
				resize_to(r);
				copy(r, *this);
			}

			return *this;
		}

		template<class RMat>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IDenseMatrix<RMat, T>& r)
		{
			resize_to(r);
			copy(r.derived(), *this);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixXpr<Expr, T>& r)
		{
			resize_to(r);
			evaluate_to(r.derived(), *this);
			return *this;
		}

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixXpr<Expr, S>& r)
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(is_implicitly_convertible<S, T>::value,
					"S is NOT implicitly-convertible to T.");
#endif
			resize_to(r);
			implicitly_cast_to(r.derived(), *this);
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
			return m_internal.ptr_col(j);
		}

		LMAT_ENSURE_INLINE pointer ptr_col(const index_type j)
		{
			return m_internal.ptr_col(j);
		}

		LMAT_ENSURE_INLINE index_type offset(const index_type i, const index_type j) const
		{
			return matrix_indexer<CTRows, CTCols>::offset(lead_dim(), i, j);
		}

		LMAT_ENSURE_INLINE const_reference elem(const index_type i, const index_type j) const
		{
			return m_internal.ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE reference elem(const index_type i, const index_type j)
		{
			return m_internal.ptr_data()[offset(i, j)];
		}

		LMAT_ENSURE_INLINE const_reference operator[] (const index_type i) const
		{
			return m_internal.ptr_data()[i];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			return m_internal.ptr_data()[i];
		}

		LMAT_ENSURE_INLINE void resize(index_type m, index_type n)
		{
			m_internal.resize(m, n);
		}

	private:
		template<class RMat>
		LMAT_ENSURE_INLINE void resize_to(const IMatrixXpr<RMat, T>& r)
		{
			if (!has_same_size(*this, r))
				resize(r.nrows(), r.ncolumns());
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

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_col(const IMatrixXpr<Expr, S>& r) : base_mat_t(r) { }

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

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_col& operator = (const IMatrixXpr<Expr, S>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE
		void resize(index_t m, index_t n)
		{
			base_mat_t::resize(m, n);
		}

		LMAT_ENSURE_INLINE
		void resize(index_t n)
		{
			base_mat_t::resize(n, 1);
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

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_row(const IMatrixXpr<Expr, S>& r) : base_mat_t(r) { }

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

		template<class Expr, typename S>
		LMAT_ENSURE_INLINE dense_row& operator = (const IMatrixXpr<Expr, S>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE
		void resize(index_t m, index_t n)
		{
			base_mat_t::resize(m, n);
		}

		LMAT_ENSURE_INLINE
		void resize(index_t n)
		{
			base_mat_t::resize(1, n);
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
