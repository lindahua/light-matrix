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

#include <light_mat/math/simd_packs.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{

	// type inference devices

	namespace meta
	{
		template<class Kernel>
		struct kernel_value_type
		{
			typedef typename Kernel::value_type type;
		};

		template<class Fun>
		struct fun_value_type
		{
			typedef typename Fun::value_type type;
		};
	}


	// access tags

	namespace atags
	{
		// access units

		struct scalar { };

		template<typename Kind>
		struct simd { };

		// access patterns

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


	typedef atags::simd<default_simd_kind> default_access_unit_t;

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


	template<class Arg, typename ATag>
	LMAT_ENSURE_INLINE
	inline in_wrap<Arg, ATag> in_(const Arg& arg, ATag)
	{
		return in_wrap<Arg, ATag>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline in_wrap<Arg, atags::normal> in_(const Arg& arg)
	{
		return in_wrap<Arg, atags::normal>(arg);
	}

	template<class Arg, typename ATag>
	LMAT_ENSURE_INLINE
	inline out_wrap<Arg, ATag> out_(Arg& arg, ATag)
	{
		return out_wrap<Arg, ATag>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline out_wrap<Arg, atags::normal> out_(Arg& arg)
	{
		return out_wrap<Arg, atags::normal>(arg);
	}

	template<class Arg, typename ATag>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<Arg, ATag> in_out_(Arg& arg, ATag)
	{
		return in_out_wrap<Arg, ATag>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline in_out_wrap<Arg, atags::normal> in_out_(Arg& arg)
	{
		return in_out_wrap<Arg, atags::normal>(arg);
	}

}

#endif 
