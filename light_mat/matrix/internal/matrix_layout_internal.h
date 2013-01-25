/**
 * @file matrix_layout_internal.h
 *
 * Internal implementation of matrix layout
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_LAYOUT_INTERNAL_H_
#define LIGHTMAT_MATRIX_LAYOUT_INTERNAL_H_

#include <light_mat/matrix/matrix_shape.h>

namespace lmat { namespace internal {

	inline index_t raise_no_linear_offset()
	{
		throw invalid_operation("Linear offset is only supported for compile-time vectors");
	}

	// contiguous layout

	template<index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline index_t cont_cm_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j)
	{
		return i + shape.nrows() * j;
	}

	template<index_t M>
	LMAT_ENSURE_INLINE
	inline index_t cont_cm_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j)
	{
		return i;
	}

	template<index_t N>
	LMAT_ENSURE_INLINE
	inline index_t cont_cm_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j)
	{
		return shape.nrows() * j;
	}

	LMAT_ENSURE_INLINE
	inline index_t cont_cm_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j)
	{
		return 0;
	}


	// block layout

	template<index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j, index_t ldim)
	{
		return i + ldim * j;
	}

	template<index_t M>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j, index_t ldim)
	{
		return i;
	}

	template<index_t N>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j, index_t ldim)
	{
		return ldim * j;
	}

	LMAT_ENSURE_INLINE
	inline index_t block_cm_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j, index_t ldim)
	{
		return 0;
	}

	template<index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_linoffset(const matrix_shape<M, N>& shape, index_t i, index_t ldim)
	{
		return raise_no_linear_offset();
	}

	template<index_t M>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_linoffset(const matrix_shape<M, 1>& shape, index_t i, index_t ldim)
	{
		return i;
	}

	template<index_t N>
	LMAT_ENSURE_INLINE
	inline index_t block_cm_linoffset(const matrix_shape<1, N>& shape, index_t i, index_t ldim)
	{
		return ldim * i;
	}

	LMAT_ENSURE_INLINE
	inline index_t block_cm_linoffset(const matrix_shape<1, 1>& shape, index_t i, index_t ldim)
	{
		return 0;
	}


	// grid layout

	template<index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline index_t grid_sub2offset(const matrix_shape<M, N>& shape, index_t i, index_t j, index_t rs, index_t cs)
	{
		return rs * i + cs * j;
	}

	template<index_t M>
	LMAT_ENSURE_INLINE
	inline index_t grid_sub2offset(const matrix_shape<M, 1>& shape, index_t i, index_t j, index_t rs, index_t cs)
	{
		return rs * i;
	}

	template<index_t N>
	LMAT_ENSURE_INLINE
	inline index_t grid_sub2offset(const matrix_shape<1, N>& shape, index_t i, index_t j, index_t rs, index_t cs)
	{
		return cs * j;
	}

	LMAT_ENSURE_INLINE
	inline index_t grid_sub2offset(const matrix_shape<1, 1>& shape, index_t i, index_t j, index_t rs, index_t cs)
	{
		return 0;
	}

	template<index_t M, index_t N>
	LMAT_ENSURE_INLINE
	inline index_t grid_linoffset(const matrix_shape<M, N>& shape, index_t i, index_t rs, index_t cs)
	{
		return raise_no_linear_offset();
	}

	template<index_t M>
	LMAT_ENSURE_INLINE
	inline index_t grid_linoffset(const matrix_shape<M, 1>& shape, index_t i, index_t rs, index_t cs)
	{
		return i * rs;
	}

	template<index_t N>
	LMAT_ENSURE_INLINE
	inline index_t grid_linoffset(const matrix_shape<1, N>& shape, index_t i, index_t rs, index_t cs)
	{
		return i * cs;
	}

	LMAT_ENSURE_INLINE
	inline index_t grid_linoffset(const matrix_shape<1, 1>& shape, index_t i, index_t rs, index_t cs)
	{
		return 0;
	}

} }

#endif 
