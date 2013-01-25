/**
 * @file macc_policy.h
 *
 * @brief Matrix access policy
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MACC_POLICY_H_
#define LIGHTMAT_MACC_POLICY_H_

#include <light_mat/mateval/mateval_fwd.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// Policies

	struct linear_ { };
	struct percol_ { };

	template<typename Acc, typename U> struct macc_ { };

	template<typename U>
	LMAT_ENSURE_INLINE
	inline bool use_linear_acc(macc_<linear_, U>)
	{
		return true;
	}

	template<typename U>
	LMAT_ENSURE_INLINE
	inline bool use_linear_acc(macc_<percol_, U>)
	{
		return false;
	}

	template<typename Acc, typename U>
	LMAT_ENSURE_INLINE
	inline bool use_simd(macc_<Acc, U>)
	{
		return false;
	}

	template<typename Acc, typename Kind>
	LMAT_ENSURE_INLINE
	inline bool use_simd(macc_<Acc, simd_<Kind> >)
	{
		return true;
	}


	/********************************************
	 *
	 *  Linear Access support
	 *
	 ********************************************/

	namespace internal
	{
		template<typename Mat>
		struct _matrix_supports_linear_macc :
		public meta::and_<
			meta::is_regular_mat<Mat>,
			meta::supports_linear_index<Mat>
		> { };
	}

	template<typename A>
	struct supports_linear_access
	: public internal::_matrix_supports_linear_macc<A> { };

	template<typename T, typename ATag>
	struct supports_linear_access<arg_wrap<T, ATag> > : public meta::false_ { };

	template<typename T>
	struct supports_linear_access<arg_wrap<T, atags::single> >
	: public meta::true_ { };

	template<typename A>
	struct supports_linear_access<arg_wrap<A, atags::in> >
	: public supports_linear_access<A> { };

	template<typename A>
	struct supports_linear_access<arg_wrap<A, atags::out> >
	: public supports_linear_access<A> { };

	template<typename A>
	struct supports_linear_access<arg_wrap<A, atags::in_out> >
	: public supports_linear_access<A> { };

	template<typename T>
	struct supports_linear_access<arg_wrap<T, atags::sum> >
	: public meta::true_ { };

	template<typename T>
	struct supports_linear_access<arg_wrap<T, atags::max> >
	: public meta::true_ { };

	template<typename T>
	struct supports_linear_access<arg_wrap<T, atags::min> >
	: public meta::true_ { };


	/********************************************
	 *
	 *  SIMD support
	 *
	 ********************************************/

	template<typename A, typename Kind> struct supports_simd;

	namespace internal
	{
		template<typename Mat, bool IsRegular, typename Kind>
		struct _matrix_supports_simd : public meta::false_ { };

		template<typename Mat, typename Kind>
		struct _matrix_supports_simd<Mat, true, Kind>
		{
			typedef typename meta::value_type_of<Mat>::type VT;

			static const bool _bt = supports_simd<VT, Kind>::value;
			static const bool _bs = meta::is_contiguous<Mat>::value ||
					(meta::is_percol_contiguous<Mat>::value && !meta::is_row<Mat>::value);

			static const bool value = _bt && _bs;
		};
	}

	template<typename A, typename Kind>
	struct supports_simd
	: public internal::_matrix_supports_simd<A, meta::is_regular_mat<A>::value, Kind> { };

	template<>
	struct supports_simd<float, sse_t> : public meta::true_ { };

	template<>
	struct supports_simd<float, avx_t> : public meta::true_ { };

	template<>
	struct supports_simd<double, sse_t> : public meta::true_ { };

	template<>
	struct supports_simd<double, avx_t> : public meta::true_ { };

	template<typename A, typename ATag, typename Kind>
	struct supports_simd<arg_wrap<A, ATag>, Kind>
	: public supports_simd<A, Kind> { };


	/********************************************
	 *
	 *  preferred policy
	 *
	 ********************************************/

	namespace internal
	{
		template<class Kernel, typename Kind, bool IsSimdizable>
		struct _kernel_packwidth
		{
			static const unsigned int value = 1;
		};

		template<class Kernel, typename Kind>
		struct _kernel_packwidth<Kernel, Kind, true>
		{
			typedef typename Kernel::value_type T;
			static const unsigned int value = simd_traits<T, Kind>::pack_width;
		};

		template<class Shape, class Kernel, typename... Args>
		struct _preferred_macc_policy_deriv
		{
			static const bool supp_linear =
					meta::all_<supports_linear_access<Args>...>::value;

			typedef default_simd_kind skind;

			static const bool ker_simdizable = is_simdizable<Kernel, skind>::value;

			static const bool args_supp_simd =
					meta::all_<supports_simd<Args, skind>...>::value;

			static const bool use_linear = supp_linear;

			static const index_t len = use_linear ?
					(Shape::ct_nrows * Shape::ct_ncols) : Shape::ct_nrows;

			static const unsigned int pack_width =
					internal::_kernel_packwidth<Kernel, skind, ker_simdizable>::value;

			static const bool use_simd =
					ker_simdizable &&
					args_supp_simd &&
					((unsigned int)len % pack_width == 0);
		};
	}



	template<class Shape, class Kernel, typename... Args>
	struct preferred_macc_policy
	{
		// derivation:

		typedef internal::_preferred_macc_policy_deriv<Shape, Kernel, Args...> _deriv;

		typedef typename _deriv::skind skind;
		static const bool use_linear = _deriv::use_linear;
		static const bool use_simd = _deriv::use_simd;

		// result:

		typedef typename std::conditional<use_linear, linear_, percol_>::type access;
		typedef typename std::conditional<use_simd, simd_<skind>, scalar_>::type unit;

		typedef macc_<access, unit> type;
	};


	template<index_t CM, index_t CN, class Kernel, typename... Args>
	typename preferred_macc_policy<matrix_shape<CM, CN>, Kernel, Args...>::type
	get_preferred_macc_policy(const matrix_shape<CM, CN>&, const Kernel&, const Args&...)
	{
		typedef typename preferred_macc_policy<matrix_shape<CM, CN>, Kernel, Args...>::type policy_t;
		return policy_t();
	}

	template<class Kernel, typename... Args>
	typename preferred_macc_policy<matrix_shape<0, 0>, Kernel, Args...>::type
	get_preferred_macc_policy(index_t, index_t, const Kernel&, const Args&...)
	{
		typedef typename preferred_macc_policy<matrix_shape<0, 0>, Kernel, Args...>::type policy_t;
		return policy_t();
	}


	template<typename T, class Expr, typename Dst>
	typename preferred_macc_policy<
		typename meta::common_shape<Expr, Dst>::type, copy_kernel<T>, Expr, Dst>::type
	get_preferred_expr_macc_policy(const IMatrixXpr<Expr, T>& expr, const IRegularMatrix<Dst, T>& dst)
	{
		const Expr& s = expr.derived();
		const Dst& d = dst.derived();

		return get_preferred_macc_policy(common_shape(s, d), copy_kernel<T>(), s, d);
	}


}

#endif /* MACC_POLICY_H_ */
