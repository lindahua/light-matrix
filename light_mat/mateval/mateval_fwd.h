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

#include <light_mat/math/simd_base.h>
#include <light_mat/matrix/matrix_concepts.h>

namespace lmat
{
	// access tags

	namespace atags
	{
		// access units

		struct scalar { };

		template<typename T, typename Kind>
		struct simd { };

		// access patterns

		struct normal { };
		struct repcol { };
		struct reprow { };
		struct single { };

		struct accum_col { };
		struct accum_row { };
		struct accum_single { };
	};

	// argument wrapper

	template<class Arg, typename ATag>
	class in_wrap
	{
	public:
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
		LMAT_ENSURE_INLINE
		out_wrap(Arg& a) : m_arg(a) { }

		LMAT_ENSURE_INLINE
		Arg& arg() const { return m_arg; }

	private:
		Arg& m_arg;
	};


	template<class Arg, typename ATag>
	LMAT_ENSURE_INLINE
	in_wrap<Arg, ATag> in_(const Arg& arg, ATag)
	{
		return in_wrap<Arg, ATag>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	in_wrap<Arg, atags::normal> in_(const Arg& arg)
	{
		return in_wrap<Arg, atags::normal>(arg, atags::normal);
	}

	template<class Arg, typename ATag>
	LMAT_ENSURE_INLINE
	out_wrap<Arg, ATag> out_(Arg& arg, ATag)
	{
		return out_wrap<Arg, ATag>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	out_wrap<Arg, atags::normal> out_(const Arg& arg)
	{
		return out_wrap<Arg, atags::normal>(arg, atags::normal);
	}


	struct copy_kernel
	{
		template<typename T>
		LMAT_ENSURE_INLINE
		void operator() (const T& s, T& d) const
		{
			d = s;
		}
	};

}

#endif 
