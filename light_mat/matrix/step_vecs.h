/**
 * @file step_vecs.h
 *
 * @brief Step vector classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_STEP_VECS_H_
#define LIGHTMAT_STEP_VECS_H_

#include <light_mat/matrix/ref_block.h>
#include <light_mat/matrix/ref_grid.h>

namespace lmat
{

	/********************************************
	 *
	 *  step column
	 *
	 ********************************************/

	template<typename T, index_t CM>
	class cstep_col: public cref_grid<T, CM, 1>
	{
		typedef cref_grid<T, CM, 1> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		cstep_col(const T* pdata, index_type m, index_type step)
		: base_mat_t(pdata, m, 1, step, 0) { }

		LMAT_ENSURE_INLINE
		cstep_col(const base_mat_t& s)
		: base_mat_t(s) { }

	};

	template<typename T, index_t CM>
	class step_col: public ref_grid<T, CM, 1>
	{
		typedef ref_grid<T, CM, 1> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		step_col(T* pdata, index_type m, index_type step)
		: base_mat_t(pdata, m, 1, step, 0) { }

		LMAT_ENSURE_INLINE
		step_col(const base_mat_t& s)
		: base_mat_t(s) { }

		template<class Expr>
		LMAT_ENSURE_INLINE step_col& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}
	};



	/********************************************
	 *
	 *  step row
	 *
	 ********************************************/

	template<typename T, index_t CN>
	class cstep_row: public cref_block<T, 1, CN>
	{
		typedef cref_block<T, 1, CN> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		cstep_row(const T* pdata, index_type n, index_type step)
		: base_mat_t(pdata, 1, n, step) { }

		LMAT_ENSURE_INLINE
		cstep_row(const base_mat_t& s)
		: base_mat_t(s) { }
	};

	template<typename T, index_t CN>
	class step_row: public ref_block<T, 1, CN>
	{
		typedef ref_block<T, 1, CN> base_mat_t;

	public:
		typedef index_t index_type;

		LMAT_ENSURE_INLINE
		step_row(T* pdata, index_type n, index_type step)
		: base_mat_t(pdata, 1, n, step) { }

		LMAT_ENSURE_INLINE
		step_row(const base_mat_t& s)
		: base_mat_t(s) { }

		LMAT_ENSURE_INLINE step_row& operator = (const base_mat_t& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}

		template<class Expr>
		LMAT_ENSURE_INLINE step_row& operator = (const IMatrixXpr<Expr, T>& r)
		{
			base_mat_t::operator = (r);
			return *this;
		}
	};

}

#endif /* STEP_VECS_H_ */
