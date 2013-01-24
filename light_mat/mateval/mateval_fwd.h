/**
 * @file mateval_fwd.h
 *
 * Forward declarations for matrix evaluation module
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATEVAL_FWD_H_
#define LIGHTMAT_MATEVAL_FWD_H_

#include <light_mat/simd/simd.h>
#include <light_mat/matrix/matrix_concepts.h>
#include <light_mat/math/functor_base.h>

#define _LMAT_DEFINE_READONLY_ARG_WRAPPER(TagName) \
	template<class Arg> \
	class arg_wrap<Arg, atags::TagName> \
	{ \
	public: \
		typedef Arg arg_type; \
		typedef atags::TagName tag_type; \
		LMAT_ENSURE_INLINE \
		arg_wrap(const Arg& a) : m_arg(a) { } \
		LMAT_ENSURE_INLINE \
		const Arg& arg() const { return m_arg; } \
	private: \
		const Arg& m_arg; \
	};

#define _LMAT_DEFINE_WRITABLE_ARG_WRAPPER(TagName) \
	template<class Arg> \
	class arg_wrap<Arg, atags::TagName> \
	{ \
	public: \
		typedef Arg arg_type; \
		typedef atags::TagName tag_type; \
		LMAT_ENSURE_INLINE \
		arg_wrap(Arg& a) : m_arg(a) { } \
		LMAT_ENSURE_INLINE \
		Arg& arg() const { return m_arg; } \
	private: \
		Arg& m_arg; \
	};


#define _LMAT_DEFINE_READONLY_ARGWRAP_FUN(TagName, WFun) \
	template<class Arg, typename T> \
	LMAT_ENSURE_INLINE \
	inline arg_wrap<Arg, atags::TagName> WFun(const IEWiseMatrix<Arg, T>& arg) \
	{ \
		return arg_wrap<Arg, atags::TagName>(arg.derived()); \
	}

#define _LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(TagName, WFun) \
	template<class Arg, typename T> \
	LMAT_ENSURE_INLINE \
	inline arg_wrap<Arg, atags::TagName> WFun(IRegularMatrix<Arg, T>& arg) \
	{ \
		return arg_wrap<Arg, atags::TagName>(arg.derived()); \
	}


namespace lmat
{
	// access units

	struct scalar_ { };

	template<typename Kind>
	struct simd_ { };

	typedef simd_<default_simd_kind> default_access_unit_t;

	// access tags

	namespace atags
	{
		struct in { };
		struct out { };
		struct in_out { };

		struct single { };
		struct repcol { };
		struct reprow { };

		struct sum { };
		struct max { };
		struct min { };

		struct colwise_sum { };
		struct colwise_max { };
		struct colwise_min { };

		struct rowwise_sum { };
		struct rowwise_max { };
		struct rowwise_min { };
	}

	// argument wrapper

	template<class Arg, typename ATag>
	class arg_wrap;

	_LMAT_DEFINE_READONLY_ARG_WRAPPER(in)
	_LMAT_DEFINE_READONLY_ARG_WRAPPER(single)
	_LMAT_DEFINE_READONLY_ARG_WRAPPER(repcol)
	_LMAT_DEFINE_READONLY_ARG_WRAPPER(reprow)

	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(out)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(in_out)

	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(sum)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(max)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(min)

	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(colwise_sum)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(colwise_max)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(colwise_min)

	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(rowwise_sum)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(rowwise_max)
	_LMAT_DEFINE_WRITABLE_ARG_WRAPPER(rowwise_min)


	// matrix-wrappers

	_LMAT_DEFINE_READONLY_ARGWRAP_FUN(in, in_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(out, out_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(in_out, in_out_)

	template<typename T>
	LMAT_ENSURE_INLINE
	inline arg_wrap<T, atags::single> const_(const T& v)
	{
		return arg_wrap<T, atags::single>(v);
	}

	_LMAT_DEFINE_READONLY_ARGWRAP_FUN(repcol, repcol_)
	_LMAT_DEFINE_READONLY_ARGWRAP_FUN(reprow, reprow_)

	template<typename T>
	LMAT_ENSURE_INLINE
	inline arg_wrap<T, atags::sum> sum_to_(T& v)
	{
		return arg_wrap<T, atags::sum>(v);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline arg_wrap<T, atags::max> max_to_(T& v)
	{
		return arg_wrap<T, atags::max>(v);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline arg_wrap<T, atags::min> min_to_(T& v)
	{
		return arg_wrap<T, atags::min>(v);
	}

	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(colwise_sum, colwise_sum_to_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(colwise_max, colwise_max_to_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(colwise_min, colwise_min_to_)

	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(rowwise_sum, rowwise_sum_to_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(rowwise_max, rowwise_max_to_)
	_LMAT_DEFINE_WRITABLE_ARGWRAP_FUN(rowwise_min, rowwise_min_to_)

}

#endif 
