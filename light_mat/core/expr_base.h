/**
 * @file expr_base.h
 *
 * The basic devices to support the implementation of expression templates
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_EXPR_BASE_H_
#define LIGHTMAT_EXPR_BASE_H_

#include <light_mat/core/basic_defs.h>

namespace lmat
{

	/********************************************
	 *
	 *  Argument forwarding and holding
	 *
	 ********************************************/

	// forward declarations

	template<class Arg> struct ref_arg_forwarder;
	template<class Arg> struct copy_arg_forwarder;
	template<class Arg> struct cache_arg_forwarder;

	template<class Arg> class ref_arg_holder;
	template<class Arg> class copy_arg_holder;
	template<class Arg> class cache_arg_holder;

	template<class Arg> class evaluated_result;

	template<class Holder>
	struct is_arg_holder
	{
		static const bool value = false;
	};

	// cross maps between forwarders and holders

	template<class Holder> struct arg_forwarder;
	template<class Forwarder> struct arg_holder;

	template<class Arg>
	struct arg_forwarder<ref_arg_holder<Arg> >
	{
		typedef ref_arg_forwarder<Arg> type;
	};

	template<class Arg>
	struct arg_forwarder<copy_arg_holder<Arg> >
	{
		typedef copy_arg_forwarder<Arg> type;
	};

	template<class Arg>
	struct arg_forwarder<cache_arg_holder<Arg> >
	{
		typedef cache_arg_forwarder<Arg> type;
	};


	template<class Arg>
	struct arg_holder<ref_arg_forwarder<Arg> >
	{
		typedef ref_arg_holder<Arg> type;
	};

	template<class Arg>
	struct arg_holder<copy_arg_forwarder<Arg> >
	{
		typedef copy_arg_holder<Arg> type;
	};

	template<class Arg>
	struct arg_holder<cache_arg_forwarder<Arg> >
	{
		typedef cache_arg_holder<Arg> type;
	};


	// Forwarders

	template<class Derived>
	struct IArgForwarder
	{
		LMAT_CRTP_REF
	};

	template<class Arg>
	struct ref_arg_forwarder : public IArgForwarder<ref_arg_forwarder<Arg> >
	{
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit ref_arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	struct copy_arg_forwarder : public IArgForwarder<copy_arg_forwarder<Arg> >
	{
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit copy_arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	struct cache_arg_forwarder : public IArgForwarder<cache_arg_forwarder<Arg> >
	{
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit cache_arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	LMAT_ENSURE_INLINE
	ref_arg_forwarder<Arg> ref_arg(const Arg& arg)
	{
		return ref_arg_forwarder<Arg>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	copy_arg_forwarder<Arg> copy_arg(const Arg& arg)
	{
		return copy_arg_forwarder<Arg>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	cache_arg_forwarder<Arg> cache_arg(const Arg& arg)
	{
		return cache_arg_forwarder<Arg>(arg);
	}


	// Holders

	template<class Arg>
	struct is_arg_holder<ref_arg_holder<Arg> >
	{
		static const bool value = true;
	};

	template<class Arg>
	struct is_arg_holder<copy_arg_holder<Arg> >
	{
		static const bool value = true;
	};

	template<class Arg>
	struct is_arg_holder<cache_arg_holder<Arg> >
	{
		static const bool value = true;
	};


	template<class Arg>
	class ref_arg_holder
	{
	public:
		typedef Arg arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit ref_arg_holder(ref_arg_forwarder<Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const arg_type& get() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	template<class Arg>
	class copy_arg_holder
	{
	public:
		typedef Arg arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit copy_arg_holder(copy_arg_forwarder<Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const arg_type& get() const
		{
			return m_arg;
		}

	private:
		Arg m_arg;
	};


	template<class Arg>
	class cache_arg_holder
	{
	public:
		typedef typename evaluated_result<Arg>::type arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit cache_arg_holder(cache_arg_forwarder<Arg> fwd)
		: m_cache(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const arg_type& get() const
		{
			return m_cache;
		}

	private:
		arg_type m_cache;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename Spec, class Arg_Holder> struct unary_expr_map;
	template<typename Spec, class Arg1_Holder, class Arg2_Holder> struct binary_expr_map;
	template<typename Spec, class Arg1_Holder, class Arg2_Holder, class Arg3_Holder> struct ternary_expr_map;

	template<typename Spec, class Arg>
	struct unary_expr_verifier
	{
		static const bool value = false;
	};

	template<typename Spec, class Arg1, class Arg2>
	struct binary_expr_verifier
	{
		static const bool value = false;
	};

	template<typename Spec, class Arg1, class Arg2, class Arg3>
	struct ternary_expr_verifier
	{
		static const bool value = false;
	};

	template<class Expr> struct expr_short_repr;

	// expression classes

	template<class Arg_Holder>
	class unary_expr
	{
	public:
		typedef typename arg_forwarder<Arg_Holder>::type arg_forwarder;
		typedef typename Arg_Holder::arg_type arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit unary_expr(arg_forwarder arg_fwd)
		: m_arg_holder(arg_fwd)
		{ }

		LMAT_ENSURE_INLINE
		const arg_type& arg() const
		{
			return m_arg_holder.get();
		}

	private:
		Arg_Holder m_arg_holder;
	};

	template<class Arg1_Holder, class Arg2_Holder>
	class binary_expr
	{
	public:
		typedef typename arg_forwarder<Arg1_Holder>::type arg1_forwarder;
		typedef typename arg_forwarder<Arg2_Holder>::type arg2_forwarder;

		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;

	public:
		LMAT_ENSURE_INLINE
		binary_expr(arg1_forwarder arg1_fwd, arg2_forwarder arg2_fwd)
		: m_arg1_holder(arg1_fwd)
		, m_arg2_holder(arg2_fwd)
		{ }

		LMAT_ENSURE_INLINE
		const arg1_type& first_arg() const
		{
			return m_arg1_holder.get();
		}

		LMAT_ENSURE_INLINE
		const arg2_type& second_arg() const
		{
			return m_arg2_holder.get();
		}

	private:
		Arg1_Holder m_arg1_holder;
		Arg2_Holder m_arg2_holder;
	};


	template<class Arg1_Holder, class Arg2_Holder, class Arg3_Holder>
	class ternary_expr
	{
	public:
		typedef typename arg_forwarder<Arg1_Holder>::type arg1_forwarder;
		typedef typename arg_forwarder<Arg2_Holder>::type arg2_forwarder;
		typedef typename arg_forwarder<Arg3_Holder>::type arg3_forwarder;

		typedef typename Arg1_Holder::arg_type arg1_type;
		typedef typename Arg2_Holder::arg_type arg2_type;
		typedef typename Arg3_Holder::arg_type arg3_type;

	public:
		LMAT_ENSURE_INLINE
		ternary_expr(arg1_forwarder arg1_fwd, arg2_forwarder arg2_fwd, arg3_forwarder arg3_fwd)
		: m_arg1_holder(arg1_fwd)
		, m_arg2_holder(arg2_fwd)
		, m_arg3_holder(arg3_fwd)
		{ }

		LMAT_ENSURE_INLINE
		const arg1_type& first_arg() const
		{
			return m_arg1_holder.get();
		}

		LMAT_ENSURE_INLINE
		const arg2_type& second_arg() const
		{
			return m_arg2_holder.get();
		}

		LMAT_ENSURE_INLINE
		const arg3_type& third_arg() const
		{
			return m_arg3_holder.get();
		}

	private:
		Arg1_Holder m_arg1_holder;
		Arg2_Holder m_arg2_holder;
		Arg3_Holder m_arg3_holder;
	};


	// functions to make expression

	template<typename Spec, class Arg_Holder>
	LMAT_ENSURE_INLINE
	typename enable_if<
		unary_expr_verifier<Spec,
			typename Arg_Holder::arg_type>,
		typename unary_expr_map<Spec,
			Arg_Holder>::type >::type
	make_expr(const Spec& spec, const typename arg_forwarder<Arg_Holder>::type& arg_fwd)
	{
		return unary_expr_map<Spec, Arg_Holder>::get(spec, arg_fwd);
	}


	template<typename Spec, class Arg1_Holder, class Arg2_Holder>
	LMAT_ENSURE_INLINE
	typename enable_if<
		binary_expr_verifier<Spec,
			typename Arg1_Holder::arg_type,
			typename Arg2_Holder::arg_type>,
		typename binary_expr_map<Spec,
			Arg1_Holder,
			Arg2_Holder>::type >::type
	make_expr(const Spec& spec,
			const typename arg_forwarder<Arg1_Holder>::type& arg1_fwd,
			const typename arg_forwarder<Arg2_Holder>::type& arg2_fwd)
	{
		return binary_expr_map<Spec, Arg1_Holder, Arg2_Holder>::get(spec, arg1_fwd, arg2_fwd);
	}


	template<typename Spec, class Arg1_Holder, class Arg2_Holder, class Arg3_Holder>
	LMAT_ENSURE_INLINE
	typename enable_if<
		ternary_expr_verifier<Spec,
			typename Arg1_Holder::arg_type,
			typename Arg2_Holder::arg_type,
			typename Arg3_Holder::arg_type>,
		typename ternary_expr_map<Spec,
			Arg1_Holder,
			Arg2_Holder,
			Arg3_Holder>::type >::type
	make_expr(const Spec& spec,
			const typename arg_forwarder<Arg1_Holder>::type& arg1_fwd,
			const typename arg_forwarder<Arg2_Holder>::type& arg2_fwd,
			const typename arg_forwarder<Arg3_Holder>::type& arg3_fwd)
	{
		return binary_expr_map<Spec, Arg1_Holder, Arg2_Holder>::get(spec, arg1_fwd, arg2_fwd, arg3_fwd);
	}


}

#endif /* EXPR_BASE_H_ */
