/*
 * @file dense_matrix.h
 *
 * The matrix to represent a dense matrix.
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_DENSE_MATRIX_H_
#define LIGHTMAT_DENSE_MATRIX_H_

#include <light_mat/matrix/regular_mat_base.h>
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
		static const int ct_num_rows = CTRows;
		static const int ct_num_cols = CTCols;

		static const bool is_readonly = false;

		typedef T value_type;
		typedef continuous_layout_cm<CTRows, CTCols> layout_type;
		typedef cpu_domain domain;
	};


	/********************************************
	 *
	 *  storage
	 *
	 ********************************************/

	namespace internal
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
			index_t capacity() const
			{
				return m_block.nelems();
			}

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
			dense_mat_storage(const dense_mat_storage& r)
			: m_block(r.m_block) { }

			LMAT_ENSURE_INLINE
			dense_mat_storage(dense_mat_storage&& r)
			: m_block(std::move(r.m_block)) { }

			LMAT_ENSURE_INLINE
			dense_mat_storage& operator = (const dense_mat_storage& r)
			{
				if (this != &r)
				{
					m_block = r.m_block;
				}
				return *this;
			}

			LMAT_ENSURE_INLINE
			dense_mat_storage& operator = (dense_mat_storage&& r)
			{
				if (this != &r)
				{
					m_block = std::move(r.m_block);
				}
				return *this;
			}

			LMAT_ENSURE_INLINE
			const T* pdata() const { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			T* pdata() { return m_block.ptr_data(); }

			LMAT_ENSURE_INLINE
			index_t capacity() const
			{
				return m_block.nelems();
			}

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
	class dense_matrix : public regular_mat_base<dense_matrix<T, CTRows, CTCols>, T>
	{
		static_assert(meta::is_supported_matrix_value_type<T>::value,
				"T must be a supported matrix value type.");

		static_assert(CTRows >= 0 && CTCols >= 0,
				"CTRows and CTCols must be non-negative numbers.");

	public:
		LMAT_DEFINE_REGMAT_TYPES(T)
		typedef continuous_layout_cm<CTRows, CTCols> layout_type;

	public:
		LMAT_ENSURE_INLINE dense_matrix()
		: m_layout()
		, m_store()
		{
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n)
		: m_layout(m, n)
		, m_store(m_layout.nelems())
		{
		}

		template<class Setter>
		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n, const IMemorySetter<Setter>& setter)
		: m_layout(m, n)
		, m_store(m_layout.nelems())
		{
			apply(setter.derived(), m * n, m_store.pdata());
		}

		LMAT_ENSURE_INLINE dense_matrix(const dense_matrix& s)
		: m_layout(s.m_layout)
		, m_store(s.m_store)
		{
		}

		LMAT_ENSURE_INLINE dense_matrix(dense_matrix&& s)
		: m_layout(std::move(s.m_layout))
		, m_store(std::move(s.m_store))
		{
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix(const IMatrixXpr<Expr, T>& r)
		: m_layout(r.nrows(), r.ncolumns())
		, m_store(m_layout.nelems())
		{
			evaluate(r.derived(), *this);
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n, const row_major_initializer<T>& in)
		: m_layout(m, n)
		, m_store(m_layout.nelems())
		{
			row_major_initialize(*this, in._list);
		}

		LMAT_ENSURE_INLINE dense_matrix(index_t m, index_t n, const col_major_initializer<T>& in)
		: m_layout(m, n)
		, m_store(m_layout.nelems())
		{
			col_major_initialize(*this, in._list);
		}

		LMAT_ENSURE_INLINE void swap(dense_matrix& s)
		{
			using std::swap;

			swap(m_layout, s.m_layout);
			m_store.swap(s.m_store);
		}

	public:

		LMAT_ENSURE_INLINE dense_matrix& operator = (const dense_matrix& r)
		{
			if (this != &r)
			{
				assign(r);
			}
			return *this;
		}

		LMAT_ENSURE_INLINE dense_matrix& operator = (dense_matrix&& r)
		{
			if (this != &r)
			{
				m_layout = std::move(r.m_layout);
				m_store = std::move(r.m_store);
			}
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE dense_matrix& operator = (const IMatrixXpr<Expr, T>& r)
		{
			assign(r);
			return *this;
		}

	public:
		LMAT_ENSURE_INLINE const layout_type& layout() const
		{
			return m_layout;
		}

		LMAT_ENSURE_INLINE const_pointer ptr_data() const
		{
			return m_store.pdata();
		}

		LMAT_ENSURE_INLINE pointer ptr_data()
		{
			return m_store.pdata();
		}

		LMAT_ENSURE_INLINE void require_size(index_t m, index_t n)
		{
			if (!(m == this->nrows() && n == this->ncolumns()))
			{
				layout_type new_layout(m, n);

				if (m_store.capacity() != new_layout.nelems())
				{
					storage_t new_store(new_layout.nelems());
					m_store.swap(new_store);
				}
				m_layout = new_layout;
			}
		}

	private:
		template<class Expr>
		LMAT_ENSURE_INLINE void assign(const IMatrixXpr<Expr, T>& r)
		{
			require_size(r.nrows(), r.ncolumns());
			evaluate(r.derived(), *this);
		}

	private:
		typedef internal::dense_mat_storage<T, CTRows * CTCols> storage_t;

		layout_type m_layout;
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

		template<class Setter>
		LMAT_ENSURE_INLINE
		dense_col(index_t m, const IMemorySetter<Setter>& setter) : base_mat_t(m, 1, setter) { }

		LMAT_ENSURE_INLINE dense_col(const base_mat_t& s) : base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_col(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

		LMAT_ENSURE_INLINE dense_col(const std::initializer_list<T>& in)
		: base_mat_t(static_cast<index_t>(in.size()), 1)
		{
			vec_initialize(*this, in);
		}

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

		LMAT_ENSURE_INLINE void require_size(index_t m)
		{
			base_mat_t::require_size(m, 1);
		}
	};


	template<typename T, int CTCols>
	class dense_row : public dense_matrix<T, 1, CTCols>
	{
		typedef dense_matrix<T, 1, CTCols> base_mat_t;

	public:
		LMAT_ENSURE_INLINE dense_row() : base_mat_t(1, CTCols) { }

		LMAT_ENSURE_INLINE explicit dense_row(index_t n) : base_mat_t(1, n) { }

		template<class Setter>
		LMAT_ENSURE_INLINE
		dense_row(index_t n, const IMemorySetter<Setter>& setter) : base_mat_t(1, n, setter) { }

		LMAT_ENSURE_INLINE dense_row(const base_mat_t& s) : base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE dense_row(const IMatrixXpr<Expr, T>& r) : base_mat_t(r) { }

		LMAT_ENSURE_INLINE dense_row(const std::initializer_list<T>& in)
		: base_mat_t(1, static_cast<index_t>(in.size()))
		{
			vec_initialize(*this, in);
		}

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

		LMAT_ENSURE_INLINE void require_size(index_t n)
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
	inline dense_matrix<T, meta::nrows<Expr>::value, meta::ncols<Expr>::value>
	eval(const IMatrixXpr<Expr, T>& expr)
	{
		return dense_matrix<T, meta::nrows<Expr>::value, meta::ncols<Expr>::value>(
				expr.derived());
	}

	template<typename T, class Expr>
	LMAT_ENSURE_INLINE
	inline T to_scalar(const IMatrixXpr<Expr, T>& expr)
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
