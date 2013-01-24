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

namespace lmat
{
	struct scalar_ { };

	template<typename Kind>
	struct simd_ { };


	// access tags

	namespace atags
	{
		struct normal { };
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


	typedef simd_<default_simd_kind> default_access_unit_t;

	// argument wrapper

	template<class Arg, typename ATag>
	class in_wrap
	{
	public:
		typedef Arg arg_type;
		typedef ATag tag_type;

		LMAT_ENSURE_INLINE
		in_wrap(const Arg& a) : m_arg(a) { }

		LMAT_ENSURE_INLINE
		const Arg& arg() const { return m_arg; }

	private:
		const Arg& m_arg;
	};


	template<class Arg, typename ATag>
	class out_wrap
	{
	public:
		typedef Arg arg_type;
		typedef ATag tag_type;

		LMAT_ENSURE_INLINE
		out_wrap(Arg& a) : m_arg(a) { }

		LMAT_ENSURE_INLINE
		Arg& arg() const { return m_arg; }

	private:
		Arg& m_arg;
	};


	template<class Arg, typename ATag>
	class in_out_wrap
	{
	public:
		typedef Arg arg_type;
		typedef ATag tag_type;

		LMAT_ENSURE_INLINE
		in_out_wrap(Arg& a) : m_arg(a) { }

		LMAT_ENSURE_INLINE
		Arg& arg() const { return m_arg; }

	private:
		Arg& m_arg;
	};

	// matrix-wrappers

	template<class Arg, typename T, typename ATag>
	LMAT_ENSURE_INLINE
	inline in_wrap<Arg, ATag> in_(const IEWiseMatrix<Arg, T>& arg, ATag)
	{
		return in_wrap<Arg, ATag>(arg.derived());
	}

	template<class Arg, typename T>
	LMAT_ENSURE_INLINE
	inline in_wrap<Arg, atags::normal> in_(const IEWiseMatrix<Arg, T>& arg)
	{
		return in_wrap<Arg, atags::normal>(arg.derived());
	}

	template<class Arg, typename T, typename ATag>
	LMAT_ENSURE_INLINE
	inline out_wrap<Arg, ATag> out_(const IEWiseMatrix<Arg, T>& arg, ATag)
	{
		return out_wrap<Arg, ATag>(arg.derived());
	}

	template<class Arg, typename T>
	LMAT_ENSURE_INLINE
	inline out_wrap<Arg, atags::normal> out_(IEWiseMatrix<Arg, T>& arg)
	{
		return out_wrap<Arg, atags::normal>(arg.derived());
	}

	template<class Arg, typename T, typename ATag>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<Arg, ATag> in_out_(IEWiseMatrix<Arg, T>& arg, ATag)
	{
		return in_out_wrap<Arg, ATag>(arg.derived());
	}

	template<class Arg, typename T>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<Arg, atags::normal> in_out_(IEWiseMatrix<Arg, T>& arg)
	{
		return in_out_wrap<Arg, atags::normal>(arg.derived());
	}

	// scalar-wrappers

	template<typename T>
	LMAT_ENSURE_INLINE
	inline in_wrap<T, atags::single> in_(const T& v, atags::single)
	{
		return in_wrap<T, atags::single>(v);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<T, atags::sum> in_out_(T& v, atags::sum)
	{
		return in_out_wrap<T, atags::sum>(v);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<T, atags::max> in_out_(T& v, atags::max)
	{
		return in_out_wrap<T, atags::max>(v);
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<T, atags::min> in_out_(T& v, atags::min)
	{
		return in_out_wrap<T, atags::min>(v);
	}

}

#endif 
