/**
 * @file args_alg.h
 *
 * @brief Simple algorithms on list of arguments
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_ARGS_ALG_H_
#define LIGHTMAT_ARGS_ALG_H_

#include <light_mat/common/prim_types.h>

namespace lmat
{
	/********************************************
	 *
	 *  pass
	 *
	 ********************************************/

	template<typename... T>
	LMAT_ENSURE_INLINE
	inline void pass(const T&... args) { }


	/********************************************
	 *
	 *  args_all
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0)
	{
		return a0;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1)
	{
		return a0 && a1;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2)
	{
		return a0 && a1 && a2;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3)
	{
		return a0 && a1 && a2 && a3;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4)
	{
		return a0 && a1 && a2 && a3 && a4;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5)
	{
		return a0 && a1 && a2 && a3 && a4 && a5;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6)
	{
		return a0 && a1 && a2 && a3 && a4 && a5 && a6;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7)
	{
		return a0 && a1 && a2 && a3 && a4 && a5 && a6 && a7;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7, bool a8)
	{
		return a0 && a1 && a2 && a3 && a4 && a5 && a6 && a7 && a8;
	}

	LMAT_ENSURE_INLINE
	inline bool args_all(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7, bool a8, bool a9)
	{
		return a0 && a1 && a2 && a3 && a4 && a5 && a6 && a7 && a8 && a9;
	}


	/********************************************
	 *
	 *  args_any
	 *
	 ********************************************/

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0)
	{
		return a0;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1)
	{
		return a0 || a1;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2)
	{
		return a0 || a1 || a2;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3)
	{
		return a0 || a1 || a2 || a3;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4)
	{
		return a0 || a1 || a2 || a3 || a4;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5)
	{
		return a0 || a1 || a2 || a3 || a4 || a5;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6)
	{
		return a0 || a1 || a2 || a3 || a4 || a5 || a6;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7)
	{
		return a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7, bool a8)
	{
		return a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7 || a8;
	}

	LMAT_ENSURE_INLINE
	inline bool args_any(bool a0, bool a1, bool a2, bool a3, bool a4, bool a5, bool a6, bool a7, bool a8, bool a9)
	{
		return a0 || a1 || a2 || a3 || a4 || a5 || a6 || a7 || a8 || a9;
	}


	/********************************************
	 *
	 *  args_equal
	 *
	 ********************************************/

	template<typename T0>
	LMAT_ENSURE_INLINE
	inline bool args_equal(const T0& )
	{
		return true;
	}

	template<typename T0, typename T1>
	LMAT_ENSURE_INLINE
	inline bool args_equal(const T0& a0, const T1& a1)
	{
		return a0 == a1;
	}

	template<typename T0, typename T1, typename T2>
	LMAT_ENSURE_INLINE
	inline bool args_equal(const T0& a0, const T1& a1, const T2& a2)
	{
		return a0 == a1 && a0 == a2;
	}

	template<typename T0, typename... TR>
	LMAT_ENSURE_INLINE
	inline bool args_equal(const T0& a0, const TR&... args)
	{
		return args_all((a0 == args)...);
	}

}

#endif /* ARGS_ALG_H_ */
