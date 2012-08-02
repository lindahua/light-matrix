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

	struct ref_arg_t { };
	struct copy_arg_t { };
	struct cache_arg_t { };

	template<class Arg> class evaluated_result;

	template<typename ArgHP, class Arg> struct arg_forwarder;
	template<typename ArgHP, class Arg> struct arg_holder;

	// Forwarders

	template<class Arg>
	struct arg_forwarder<ref_arg_t, Arg>
	{
		typedef Arg arg_type;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	struct arg_forwarder<copy_arg_t, Arg>
	{
		typedef Arg arg_type;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	struct arg_forwarder<cache_arg_t, Arg>
	{
		typedef Arg arg_type;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(const Arg& a) : arg(a) { }
	};

	template<class Arg>
	LMAT_ENSURE_INLINE
	arg_forwarder<ref_arg_t, Arg> ref_arg(const Arg& arg)
	{
		return arg_forwarder<ref_arg_t, Arg>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	arg_forwarder<copy_arg_t, Arg> copy_arg(const Arg& arg)
	{
		return arg_forwarder<copy_arg_t, Arg>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	arg_forwarder<cache_arg_t, Arg> cache_arg(const Arg& arg)
	{
		return arg_forwarder<cache_arg_t, Arg>(arg);
	}


	// Holders

	template<class Arg>
	class arg_holder<ref_arg_t, Arg>
	{
	public:
		typedef Arg internal_arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit arg_holder(arg_forwarder<ref_arg_t, Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const internal_arg_type& get() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	template<class Arg>
	class arg_holder<copy_arg_t, Arg>
	{
	public:
		typedef Arg internal_arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit arg_holder(arg_forwarder<copy_arg_t, Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const internal_arg_type& get() const
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
		typedef typename evaluated_result<Arg>::type internal_arg_type;

	public:
		LMAT_ENSURE_INLINE
		explicit cache_arg_holder(arg_forwarder<cache_arg_t, Arg> fwd)
		: m_cache(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const internal_arg_type& get() const
		{
			return m_cache;
		}

	private:
		internal_arg_type m_cache;
	};


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename Spec,
		typename Arg_HP, class Arg>
	struct unary_expr_map;

	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2>
	struct binary_expr_map;

	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2,
		typename Arg3_HP, class Arg3>
	struct ternary_expr_map;

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

	template<typename Arg_HP, class Arg>
	class unary_expr
	{
	public:
		typedef arg_holder<Arg_HP, Arg> arg_holder_t;
		typedef typename arg_holder_t::internal_arg_type arg_type;

		LMAT_ENSURE_INLINE
		explicit unary_expr(const arg_forwarder<Arg_HP, Arg>& arg_fwd)
		: m_arg_holder(arg_fwd)
		{ }

		LMAT_ENSURE_INLINE
		const arg_type& arg() const
		{
			return m_arg_holder.get();
		}

	private:
		arg_holder_t m_arg_holder;
	};

	template<
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2>
	class binary_expr
	{
	public:
		typedef arg_holder<Arg1_HP, Arg1> arg1_holder_t;
		typedef arg_holder<Arg2_HP, Arg2> arg2_holder_t;

		typedef typename arg1_holder_t::internal_arg_type arg1_type;
		typedef typename arg2_holder_t::internal_arg_type arg2_type;

	public:
		LMAT_ENSURE_INLINE
		binary_expr(
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
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
		arg1_holder_t m_arg1_holder;
		arg2_holder_t m_arg2_holder;
	};


	template<
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2,
		typename Arg3_HP, class Arg3>
	class ternary_expr
	{
	public:
		typedef arg_holder<Arg1_HP, Arg1> arg1_holder_t;
		typedef arg_holder<Arg2_HP, Arg2> arg2_holder_t;
		typedef arg_holder<Arg3_HP, Arg3> arg3_holder_t;

		typedef typename arg1_holder_t::internal_arg_type arg1_type;
		typedef typename arg2_holder_t::internal_arg_type arg2_type;
		typedef typename arg3_holder_t::internal_arg_type arg3_type;

	public:
		LMAT_ENSURE_INLINE
		ternary_expr(
				const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
				const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd,
				const arg_forwarder<Arg3_HP, Arg3>& arg3_fwd)
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
		arg1_holder_t m_arg1_holder;
		arg2_holder_t m_arg2_holder;
		arg3_holder_t m_arg3_holder;
	};


	// functions to make expression

	template<typename Spec,
		typename Arg_HP, class Arg>
	LMAT_ENSURE_INLINE
	typename enable_if<
		unary_expr_verifier<Spec, Arg>,
		typename unary_expr_map<Spec,
			Arg_HP, Arg
		>::type
	>::type
	make_expr(const Spec& spec,
			const arg_forwarder<Arg_HP, Arg>& arg_fwd)
	{
		return unary_expr_map<Spec,
				Arg_HP, Arg>::get(spec, arg_fwd);
	}


	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2>
	LMAT_ENSURE_INLINE
	typename enable_if<
		binary_expr_verifier<Spec, Arg1, Arg2>,
		typename binary_expr_map<Spec,
			Arg1_HP, Arg1,
			Arg2_HP, Arg2
		>::type
	>::type
	make_expr(const Spec& spec,
			const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
			const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd)
	{
		return binary_expr_map<Spec,
				Arg1_HP, Arg1,
				Arg2_HP, Arg2>::get(spec, arg1_fwd, arg2_fwd);
	}


	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2,
		typename Arg3_HP, class Arg3>
	LMAT_ENSURE_INLINE
	typename enable_if<
		ternary_expr_verifier<Spec, Arg1, Arg2, Arg3>,
		typename ternary_expr_map<Spec,
			Arg1_HP, Arg1,
			Arg2_HP, Arg2,
			Arg3_HP, Arg3
		>::type
	>::type
	make_expr(const Spec& spec,
			const arg_forwarder<Arg1_HP, Arg1>& arg1_fwd,
			const arg_forwarder<Arg2_HP, Arg2>& arg2_fwd,
			const arg_forwarder<Arg3_HP, Arg3>& arg3_fwd)
	{
		return ternary_expr_map<Spec,
				Arg1_HP, Arg1,
				Arg2_HP, Arg2,
				Arg3_HP, Arg3>::get(spec, arg1_fwd, arg2_fwd, arg3_fwd);
	}


}

#endif /* EXPR_BASE_H_ */
