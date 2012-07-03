/*
 * @file functor_base.h
 *
 * The basic definitions for functors
 *
 * @author Dahua Lin
 */

#ifndef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_FUNCTOR_BASE_H_
#define LIGHTMAT_FUNCTOR_BASE_H_

#include <light_mat/core/basic_defs.h>

namespace lmat
{
	template<class Fun>
	struct is_unary_ewise_functor
	{
		static const bool value = false;
	};

	template<class Fun>
	struct is_binary_ewise_functor
	{
		static const bool value = false;
	};

	template<class Fun>
	struct supports_simd
	{
		static const bool value = false;
	};

	template<typename Arg, typename Result>
	struct unary_ewise_functor
	{
		typedef Arg arg_type;
		typedef Result result_type;
	};

	template<typename Arg1, typename Arg2, typename Result>
	struct binary_ewise_functor
	{
		typedef Arg1 first_arg_type;
		typedef Arg2 second_arg_type;
		typedef Result result_type;
	};

	template<typename T>
	struct unary_numeric_ewise_functor : public unary_ewise_functor<T, T> { };

	template<typename T>
	struct binary_numeric_ewise_functor : public binary_ewise_functor<T, T, T> { };

	template<typename T>
	struct unary_mpred_ewise_functor : public unary_ewise_functor<T, mask_t<T> > { };

	template<typename T>
	struct binary_mpred_ewise_functor : public binary_ewise_functor<T, T, mask_t<T> > { };

	template<typename T>
	struct unary_mask_ewise_functor : public unary_ewise_functor<mask_t<T>, mask_t<T> > { };

	template<typename T>
	struct binary_mask_ewise_functor : public binary_ewise_functor<mask_t<T>, mask_t<T>, mask_t<T> > { };


	// conversion

	template<typename S, typename T>
	struct type_converter : public unary_ewise_functor<S, T>
	{
		LMAT_ENSURE_INLINE
		T operator() (const S& s) const
		{
			return T(s);
		}
	};

}

// Useful macros

#define LMAT_DECLARE_AS_UNARY_EWISE_FUNCTOR( Fun, SuppSIMD ) \
	template<> \
	struct is_unary_ewise_functor< Fun > { static const bool value = true; }; \
	template<> \
	struct supports_simd< Fun > { static const bool value = SuppSIMD; };

#define LMAT_DECLARE_AS_BINARY_EWISE_FUNCTOR( Fun, SuppSIMD ) \
	template<> \
	struct is_binary_ewise_functor< Fun > { static const bool value = true; }; \
	template<> \
	struct supports_simd< Fun > { static const bool value = SuppSIMD; };


#define LMAT_DECLARE_AS_UNARY_EWISE_TFUNCTOR( FunT, SuppSIMD ) \
	template<typename T> \
	struct is_unary_ewise_functor< FunT<T> > { static const bool value = true; }; \
	template<typename T> \
	struct supports_simd< FunT<T> > { static const bool value = SuppSIMD; };

#define LMAT_DECLARE_AS_BINARY_EWISE_TFUNCTOR( FunT, SuppSIMD ) \
	template<typename T> \
	struct is_binary_ewise_functor< FunT<T> > { static const bool value = true; }; \
	template<typename T> \
	struct supports_simd< FunT<T> > { static const bool value = SuppSIMD; };


#define LMAT_DEFINE_UNARY_NUMERIC_EWISE_TFUNCTOR( FunT, ImplFun, SuppSIMD ) \
	template<typename T> \
	struct FunT : public unary_numeric_ewise_functor<T> { \
		T operator()(const T& x) const { return ImplFun(x); } \
	}; \
	template<typename T> \
	struct is_unary_ewise_functor< FunT<T> > { static const bool value = true; }; \
	template<typename T> \
	struct supports_simd< FunT<T> > { static const bool value = SuppSIMD; };

#define LMAT_DEFINE_BINARY_NUMERIC_EWISE_TFUNCTOR( FunT, ImplFun, SuppSIMD ) \
	template<typename T> \
	struct FunT : public binary_numeric_ewise_functor<T> { \
		T operator()(const T& x, const T& y) const { return ImplFun(x, y); } \
	}; \
	template<typename T> \
	struct is_binary_ewise_functor< FunT<T> > { static const bool value = true; }; \
	template<typename T> \
	struct supports_simd< FunT<T> > { static const bool value = SuppSIMD; };

#endif 
