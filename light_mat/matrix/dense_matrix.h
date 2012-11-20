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

#include <light_mat/matrix/matrix_base.h>
#include <light_mat/common/block.h>

#include <algorithm> // for std::swap

namespace lmat
{

	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols>
	struct matrix_traits<dense_matrix<T, CTRows, CTCols> >
	{
		static const int num_dimensions = 2;
		static const int compile_time_num_rows = CTRows;
		static const int compile_time_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef cpu_domain domain;
	};

	template<typename T, int CTRows, int CTCols>
	struct ct_has_continuous_layout<dense_matrix<T, CTRows, CTCols> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols>
	struct is_linear_accessible<dense_matrix<T, CTRows, CTCols> >
	{
		static const bool value = true;
	};

	template<typename T, int CTRows, int CTCols, class DMat>
	struct default_matrix_eval_policy<dense_matrix<T, CTRows, CTCols>, DMat>
	{
		typedef matrix_copy_policy type;
	};


	/********************************************
	 *
	 *  storage
	 *
	 ********************************************/

	namespace detail
	{
		template<typename T, int CTSize>
		class dense_mat_storage
		{
#ifdef LMAT_USE_STATIC_ASSERT
			static_assert(CTSize > 0, "CTSize must be a positive integer");
#endif
		public:
			LMAT_ENSURE_INLINE
			dense_mat_storage() : m_block() { }

			LMAT_ENSURE_INLINE
			dense_mat_storage(index_t siz) : m_block() { }

			LMAT_ENSURE_INLINE
			const T* pdata() const { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			T* pdata() { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			void swap(dense_mat_storage& other)
			{
				m_block.swap(other.m_block);
			}

		private:
			sblock<T, CTSize> m_block;
		};


		template<typename T>
		class dense_mat_storage<T, 0>
		{
		public:
			LMAT_ENSURE_INLINE
			dense_mat_storage() : m_block() { }

			LMAT_ENSURE_INLINE
			dense_mat_storage(index_t siz) : m_block(siz) { }

			LMAT_ENSURE_INLINE
			const T* pdata() const { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			T* pdata() { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			void swap(dense_mat_storage& other)
			{
				m_block.swap(other.m_block);
			}

		private:
			dblock<T> m_block;
		};
	}

	/********************************************
	 *
	 *  dense_matrix
	 *
	 ********************************************/

	template<typename T, int CTRows, int CTCols>
	class dense_matrix : public IDenseMatrix<dense_matrix<T, CTRows, CTCols>, T>
	{
#ifdef LMAT_USE_STATIC_ASSERT
		static_assert(is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type.");

		static_assert(CTRows >= 0 && CTCols >= 0,
				"CTRows and CTCols must be non-negative numbers.");
#endif

	public:
		LMAT_MAT_TRAITS_DEFS(T)

	public:
		LMAT_ENSURE_INLINE dense_matrix()
		: m_shape()
		, m_store()
		{
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n)
		: m_shape(m, n)
		, m_store(m_shape.nelems())
		{
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n, zero_t)
		: m_shape(m, n)
		, m_store(m_shape.nelems())
		{
			zero_mem(m * n, m_store.pdata());
		}

		template<class Setter>
		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n,
				const IMemorySetter<Setter, T>& setter)
		: m_shape(m, n)
		, m_store(m_shape.nelems())
		{
			setter.set(m * n, m_store.pdata());
		}

		LMAT_ENSURE_INLINE dense_matrix(const dense_matrix& s)
		: m_shape(s.m_shape)
		, m_store(s.m_store)
		{
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix(const IMatrixXpr<Expr, T>& r)
		: m_shape(r.nrows(), r.ncolumns())
		, m_store(m_shape.nelems())
		{
			default_evaluate(r, *this);
		}

		LMAT_ENSURE_INLINE void swap(dense_matrix& s)
		{
			using std::swap;

			swap(m_shape, s.m_shape);
			m_store.swap(s.m_store);
		}

	public:

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
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE size_type size() const
		{
			return static_cast<size_type>(nelems());
		}

		LMAT_ENSURE_INLINE index_type nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_type ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_type lead_dim() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_store.pdata();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_store.pdata();
		}

		LMAT_ENSURE_INLINE index_type offset(const index_type i, const index_type j) const
		{
			return sub2offset(column_major_layout(), m_shape, i, j);
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
			LMAT_CHECK_IDX(i, nelems())
			return ptr_data()[i];
		}

		LMAT_ENSURE_INLINE reference operator[] (const index_type i)
		{
			LMAT_CHECK_IDX(i, nelems())
			return ptr_data()[i];
		}

		LMAT_ENSURE_INLINE void require_size(index_type m, index_type n)
		{
			if (!(m == nrows() && n == ncolumns()))
			{
				matrix_shape<CTRows, CTCols> new_shape(m, n);

				if (new_shape.nelems() != nelems())
				{
					storage_t new_store(new_shape.nelems());
					m_store.swap(new_store);
				}
				m_shape = new_shape;
			}
		}

	private:
		typedef detail::dense_mat_storage<T, CTRows * CTCols> storage_t;

		matrix_shape<CTRows, CTCols> m_shape;
		storage_t m_store;
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

	template<typename T, int CTRows>
	class dense_col : public dense_matrix<T, CTRows, 1>
	{
		typedef dense_matrix<T, CTRows, 1> base_mat_t;

	public:
		LMAT_ENSURE_INLINE dense_col() : base_mat_t(CTRows, 1) { }

		LMAT_ENSURE_INLINE explicit dense_col(index_t m) : base_mat_t(m, 1) { }

		LMAT_ENSURE_INLINE dense_col(index_t m, zero_t) : base_mat_t(m, 1, zero_t()) { }

		template<class Setter>
		LMAT_ENSURE_INLINE
		dense_col(index_t m, const IMemorySetter<Setter, T>& setter) : base_mat_t(m, 1, setter) { }

		LMAT_ENSURE_INLINE dense_col(const base_mat_t& s) : base_mat_t(s) { }

		LMAT_ENSURE_INLINE dense_col(const dense_col& s) : base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_col(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

	public:

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


	template<typename T, int CTCols>
	class dense_row : public dense_matrix<T, 1, CTCols>
	{
		typedef dense_matrix<T, 1, CTCols> base_mat_t;

	public:
		LMAT_ENSURE_INLINE dense_row() : base_mat_t(1, CTCols) { }

		LMAT_ENSURE_INLINE explicit dense_row(index_t n) : base_mat_t(1, n) { }

		LMAT_ENSURE_INLINE dense_row(index_t n, zero_t) : base_mat_t(1, n, zero_t()) { }

		template<class Setter>
		LMAT_ENSURE_INLINE
		dense_row(index_t n, const IMemorySetter<Setter,T>& setter) : base_mat_t(1, n, setter) { }

		LMAT_ENSURE_INLINE dense_row(const base_mat_t& s) : base_mat_t(s) { }

		LMAT_ENSURE_INLINE dense_row(const dense_row& s) : base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_row(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

	public:

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

	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat, 0, 0)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat22, 2, 2)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat23, 2, 3)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat32, 3, 2)
	LMAT_MATRIX_TYPEDEFS2(dense_matrix, mat33, 3, 3)

	LMAT_MATRIX_TYPEDEFS1(dense_col, col, 0)
	LMAT_MATRIX_TYPEDEFS1(dense_col, col2, 2)
	LMAT_MATRIX_TYPEDEFS1(dense_col, col3, 3)

	LMAT_MATRIX_TYPEDEFS1(dense_row, row, 0)
	LMAT_MATRIX_TYPEDEFS1(dense_row, row2, 2)
	LMAT_MATRIX_TYPEDEFS1(dense_row, row3, 3)

}


#endif 
