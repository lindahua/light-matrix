/**
 * @file mat_zip.h
 *
 * @brief Matrix zip and unzip
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAT_ZIP_H_
#define LIGHTMAT_MAT_ZIP_H_

#include <light_mat/matexpr/map_expr.h>
#include <tuple>

namespace lmat
{
	/********************************************
	 *
	 *  Zip
	 *
	 ********************************************/

	struct zip_pair_ { };

	template<typename T1, typename T2>
	struct zip_fun
	{
		typedef std::pair<T1, T2> result_type;

		LMAT_ENSURE_INLINE
		result_type operator() (const T1& e1, const T2& e2) const
		{
			return std::make_pair(e1, e2);
		}
	};

	template<typename T1, typename T2>
	struct fun_map<zip_pair_, T1, T2>
	{
		typedef zip_fun<T1, T2> type;
	};

	template<typename A1, typename T1, typename A2, typename T2>
	LMAT_ENSURE_INLINE
	inline map_expr<zip_pair_, A1, A2>
	zip_pair(const IEWiseMatrix<A1, T1>& a1, const IEWiseMatrix<A2, T2>& a2)
	{
		return make_map_expr(zip_pair_(), a1, a2);
	}

	template<typename T1, typename A2, typename T2>
	LMAT_ENSURE_INLINE
	inline map_expr<zip_pair_, T1, A2>
	zip_pair_fix1(const T1& a1, const IEWiseMatrix<A2, T2>& a2)
	{
		return make_map_expr_fix1(zip_pair_(), a1, a2);
	}

	template<typename A1, typename T1, typename T2>
	LMAT_ENSURE_INLINE
	inline map_expr<zip_pair_, A1, T2>
	zip_pair_fix2(const IEWiseMatrix<A1, T1>& a1, const T2& a2)
	{
		return make_map_expr_fix2(zip_pair_(), a1, a2);
	}


	/********************************************
	 *
	 *  Zip extraction
	 *
	 ********************************************/

	template<int I> struct zip_e_ { };

	template<int I, typename T>
	struct zip_e_fun
	{
		typedef typename std::tuple_element<(size_t)I, T>::type result_type;

		LMAT_ENSURE_INLINE
		result_type operator() (const T& t) const
		{
			return std::get<(size_t)I>(t);
		}
	};

	template<int I, typename T>
	struct fun_map<zip_e_<I>, T>
	{
		typedef zip_e_fun<I, T> type;
	};

	template<typename A, typename T, int I>
	LMAT_ENSURE_INLINE
	inline map_expr<zip_e_<I>, A>
	zip_e(const IEWiseMatrix<A, T>& a, int_<I>)
	{
		return make_map_expr(zip_e_<I>(), a);
	}


	/********************************************
	 *
	 *  Unzip
	 *
	 ********************************************/

	template<typename T>
	struct unzip_kernel
	{
		template<typename E0, typename E1>
		LMAT_ENSURE_INLINE
		void operator() (const T& t, E0& e0, E1& e1) const
		{
			e0 = std::get<0>(t);
			e1 = std::get<1>(t);
		}

		template<typename E0, typename E1, typename E2>
		LMAT_ENSURE_INLINE
		void operator() (const T& t, E0& e0, E1& e1, E2& e2) const
		{
			e0 = std::get<0>(t);
			e1 = std::get<1>(t);
			e2 = std::get<2>(t);
		}
	};


	template<class S, typename T, class D0, typename E0, class D1, typename E1>
	inline void unzip(const IEWiseMatrix<S, T>& t,
			IRegularMatrix<D0, E0>& e0, IRegularMatrix<D1, E1>& e1)
	{
		const index_t m = t.nrows();
		const index_t n = t.ncolumns();

		e0.require_size(m, n);
		e1.require_size(m, n);

		ewise(unzip_kernel<T>())(t.shape(), in_(t), out_(e0), out_(e1));
	}

	template<class S, typename T, class D0, typename E0, class D1, typename E1, class D2, typename E2>
	inline void unzip(const IEWiseMatrix<S, T>& t,
			IRegularMatrix<D0, E0>& e0, IRegularMatrix<D1, E1>& e1, IRegularMatrix<D2, E2>& e2)
	{
		const index_t m = t.nrows();
		const index_t n = t.ncolumns();

		e0.require_size(m, n);
		e1.require_size(m, n);
		e2.require_size(m, n);

		ewise(unzip_kernel<T>())(t.shape(), in_(t), out_(e0), out_(e1), out_(e2));
	}


}

#endif /* ZIP_EXPR_H_ */
