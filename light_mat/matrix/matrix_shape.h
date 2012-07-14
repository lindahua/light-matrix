/**
 * @file matrix_shape.h
 *
 * The class that represents the shape of a matrix
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_SHAPE_H_
#define LIGHTMAT_MATRIX_SHAPE_H_

#include <light_mat/matrix/matrix_fwd.h>

namespace lmat
{

	/********************************************
	 *
	 *  single_dim
	 *
	 ********************************************/

	template<int N>
	struct single_dim
	{
		LMAT_ENSURE_INLINE
		single_dim() { }

		LMAT_ENSURE_INLINE
		single_dim(const index_t n)
		{
			check_arg(n == N, "single_dim: the input dimension is invalid.");
		}

		LMAT_ENSURE_INLINE
		index_t dim() const { return N; }
	};


	template<>
	struct single_dim<DynamicDim>
	{
		LMAT_ENSURE_INLINE
		single_dim() : m_dim(0) { }

		LMAT_ENSURE_INLINE
		single_dim(const index_t n) : m_dim(n)
		{ }

		LMAT_ENSURE_INLINE
		index_t dim() const { return m_dim; }

		const index_t m_dim;
	};


	/********************************************
	 *
	 *  matrix_shape
	 *
	 ********************************************/

	namespace detail
	{

		template<bool IsRow, bool IsCol> struct matrix_index_helper;

		template<> struct matrix_index_helper<false, false>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return i + j * ldim;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return i + j * ldim;
			}
		};

		template<> struct matrix_index_helper<false, true>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return i;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return i;
			}
		};

		template<> struct matrix_index_helper<true, false>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return j * ldim;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return j;
			}
		};

		template<> struct matrix_index_helper<true, true>
		{
			LMAT_ENSURE_INLINE
			static index_t offset(const index_t ldim, const index_t i, const index_t j)
			{
				return 0;
			}

			LMAT_ENSURE_INLINE
			static index_t offset_c(const index_t ldim, const index_t i, const index_t j)
			{
				return 0;
			}
		};


		template<int M, int N, bool MFixed, bool NFixed>
		struct matrix_shape_internal;

		template<>
		struct matrix_shape_internal<DynamicDim, DynamicDim, false, false>
		{
			const index_t _nrows;
			const index_t _ncols;

			LMAT_ENSURE_INLINE
			matrix_shape_internal() : _nrows(0), _ncols(0) { }

			LMAT_ENSURE_INLINE
			matrix_shape_internal(const index_t m, const index_t n)
			: _nrows(m), _ncols(n)
			{
			}

			LMAT_ENSURE_INLINE index_t nrows() const { return _nrows; }
			LMAT_ENSURE_INLINE index_t ncolumns() const { return _ncols; }
			LMAT_ENSURE_INLINE index_t nelems() const { return _nrows * _ncols; }

			LMAT_ENSURE_INLINE index_t is_empty() const { return _nrows == 0 || _ncols == 0; }
			LMAT_ENSURE_INLINE index_t is_scalar() const { return _nrows == 1 && _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_column() const { return _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_row() const { return _nrows == 1; }
			LMAT_ENSURE_INLINE index_t is_vector() const { return _nrows == 1 || _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_square() const { return _nrows == _ncols; }
		};

		template<int N>
		struct matrix_shape_internal<DynamicDim, N, false, true>
		{
			const index_t _nrows;

			LMAT_ENSURE_INLINE
			matrix_shape_internal() : _nrows(0) { }

			LMAT_ENSURE_INLINE
			matrix_shape_internal(const index_t m, const index_t n)
			: _nrows(m)
			{
				check_arg(n == N,
						"Attempted to construct a matrix with incorrect #columns.");
			}

			LMAT_ENSURE_INLINE index_t nrows() const { return _nrows; }
			LMAT_ENSURE_INLINE index_t ncolumns() const { return N; }
			LMAT_ENSURE_INLINE index_t nelems() const { return _nrows * N; }

			LMAT_ENSURE_INLINE index_t is_empty() const { return _nrows == 0; }
			LMAT_ENSURE_INLINE index_t is_scalar() const { return N == 1 && _nrows == 1; }
			LMAT_ENSURE_INLINE index_t is_column() const { return N == 1; }
			LMAT_ENSURE_INLINE index_t is_row() const { return _nrows == 1; }
			LMAT_ENSURE_INLINE index_t is_vector() const { return N == 1 || _nrows == 1; }
			LMAT_ENSURE_INLINE index_t is_square() const { return _nrows == N; }
		};

		template<int M>
		struct matrix_shape_internal<M, DynamicDim, true, false>
		{
			const index_t _ncols;

			LMAT_ENSURE_INLINE
			matrix_shape_internal() : _ncols(0) { }

			LMAT_ENSURE_INLINE
			matrix_shape_internal(const index_t m, const index_t n)
			: _ncols(n)
			{
				check_arg(m == M,
						"Attempted to construct a matrix with incorrect #rows.");
			}

			LMAT_ENSURE_INLINE index_t nrows() const { return M; }
			LMAT_ENSURE_INLINE index_t ncolumns() const { return _ncols; }
			LMAT_ENSURE_INLINE index_t nelems() const { return M * _ncols; }

			LMAT_ENSURE_INLINE index_t is_empty() const { return _ncols == 0; }
			LMAT_ENSURE_INLINE index_t is_scalar() const { return M == 1 && _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_column() const { return _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_row() const { return M == 1; }
			LMAT_ENSURE_INLINE index_t is_vector() const { return M == 1 || _ncols == 1; }
			LMAT_ENSURE_INLINE index_t is_square() const { return M == _ncols; }
		};

		template<int M, int N>
		struct matrix_shape_internal<M, N, true, true>
		{
			LMAT_ENSURE_INLINE
			matrix_shape_internal() { }

			LMAT_ENSURE_INLINE
			matrix_shape_internal(const index_t m, const index_t n)
			{
				check_arg(m == M && n == N,
						"Attempted to construct a matrix with incorrect size.");
			}

			LMAT_ENSURE_INLINE index_t nrows() const { return M; }
			LMAT_ENSURE_INLINE index_t ncolumns() const { return N; }
			LMAT_ENSURE_INLINE index_t nelems() const { return M * N; }

			LMAT_ENSURE_INLINE index_t is_empty() const { return false; }
			LMAT_ENSURE_INLINE index_t is_scalar() const { return M == 1 && N == 1; }
			LMAT_ENSURE_INLINE index_t is_column() const { return N == 1; }
			LMAT_ENSURE_INLINE index_t is_row() const { return M == 1; }
			LMAT_ENSURE_INLINE index_t is_vector() const { return M == 1 || N == 1; }
			LMAT_ENSURE_INLINE index_t is_square() const { return M == N; }
		};
	}

	template<int M, int N>
	class matrix_shape
	{
	public:
		LMAT_ENSURE_INLINE
		matrix_shape() { };

		LMAT_ENSURE_INLINE
		matrix_shape(const index_t m, const index_t n)
		: m_internal(m, n) { }

		LMAT_ENSURE_INLINE index_t nrows() const { return m_internal.nrows(); }
		LMAT_ENSURE_INLINE index_t ncolumns() const { return m_internal.ncolumns(); }
		LMAT_ENSURE_INLINE index_t nelems() const { return m_internal.nelems(); }

		LMAT_ENSURE_INLINE index_t is_empty() const { return m_internal.is_empty(); }
		LMAT_ENSURE_INLINE index_t is_scalar() const { return m_internal.is_scalar(); }
		LMAT_ENSURE_INLINE index_t is_column() const { return m_internal.is_column(); }
		LMAT_ENSURE_INLINE index_t is_row() const { return m_internal.is_row(); }
		LMAT_ENSURE_INLINE index_t is_vector() const { return m_internal.is_vector(); }
		LMAT_ENSURE_INLINE index_t is_square() const { return m_internal.is_square(); }

		LMAT_ENSURE_INLINE index_t offset(const index_t ldim, const index_t i, const index_t j)
		{
			return detail::matrix_index_helper<M==1, N==1>::offset(ldim, i, j);
		}

		LMAT_ENSURE_INLINE index_t offset_c(const index_t i, const index_t j)
		{
			return detail::matrix_index_helper<M==1, N==1>::offset(nrows(), i, j);
		}

	private:
		detail::matrix_shape_internal<M, N, (M>0), (N>0)> m_internal;
	};


	/********************************************
	 *
	 *  indexing
	 *
	 ********************************************/

}

#endif /* MATRIX_SHAPE_H_ */
