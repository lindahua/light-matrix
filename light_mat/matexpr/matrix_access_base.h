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

#include <light_mat/matexpr/matexpr_fwd.h>

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


	template<class PerColAcc, typename Ker, typename T> class percol_to_linear_accessor;

	template<class PerColAcc, typename T>
	class percol_to_linear_accessor<PerColAcc, scalar_ker, T>
			: public ILinearMatrixScalarAccessor<percol_to_linear_accessor<PerColAcc, scalar_ker, T>, T>
	{
	public:
		typedef typename percol_macc_state_map<PerColAcc>::type col_state_t;

		LMAT_ENSURE_INLINE
		percol_to_linear_accessor(const PerColAcc& acc, const col_state_t& s)
		: m_acc(acc), m_col_state(s)
		{ }


		LMAT_ENSURE_INLINE
		T get_scalar(const index_t i) const
		{
			return m_acc.get_scalar(i, m_col_state);
		}

	private:
		const PerColAcc& m_acc;
		col_state_t m_col_state;
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
	struct macc_accessor_map<scalar_expr<T>, AccCate, scalar_ker>
	{
		typedef const_scalar_accessor<T> type;
	};

}

#endif /* MATRIX_ACCESS_BASE_H_ */





