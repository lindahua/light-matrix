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

#include <light_mat/common/basic_defs.h>

// useful macros

#define LMAT_QARG( Arg_HP, Arg ) typename QArg_<Arg_HP, Arg>::type

namespace lmat
{

	/********************************************
	 *
	 *  Argument forwarding and holding
	 *
	 ********************************************/

	// forward declarations

	struct cref_arg_t { };
	struct ref_arg_t { };
	struct cpy_arg_t { };

	template<typename ArgHP, class Arg> struct arg_forwarder;
	template<typename ArgHP, class Arg> class arg_holder;

	// qualified argument

	template<class Arg> class CRefArg
	{
		typedef Arg argument;
		typedef const argument& reference;
		typedef cref_arg_t policy;
		typedef arg_forwarder<policy, Arg> forwarder;
		typedef arg_holder<policy, Arg> holder;
	};

	template<class Arg> class RefArg
	{
		typedef Arg argument;
		typedef argument& reference;
		typedef ref_arg_t policy;
		typedef arg_forwarder<policy, Arg> forwarder;
		typedef arg_holder<policy, Arg> holder;
	};

	template<class Arg> class CpyArg
	{
		typedef Arg argument;
		typedef const argument& reference;
		typedef cpy_arg_t policy;
		typedef arg_forwarder<policy, Arg> forwarder;
		typedef arg_holder<policy, Arg> holder;
	};

	template<typename Arg_HP, class Arg> struct QArg_;

	template<class Arg> struct QArg_<cref_arg_t, Arg> { typedef CRefArg<Arg> type; };
	template<class Arg> struct QArg_< ref_arg_t, Arg> { typedef  RefArg<Arg> type; };
	template<class Arg> struct QArg_< cpy_arg_t, Arg> { typedef  CpyArg<Arg> type; };


	// Forwarders

	template<class Arg>
	struct arg_forwarder<cref_arg_t, Arg>
	{
		typedef Arg arg_type;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(const Arg& a) : arg(a) { }
	};


	template<class Arg>
	struct arg_forwarder<ref_arg_t, Arg>
	{
		typedef Arg arg_type;
		Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(Arg& a) : arg(a) { }
	};


	template<class Arg>
	struct arg_forwarder<cpy_arg_t, Arg>
	{
		typedef Arg arg_type;
		const Arg& arg;

		LMAT_ENSURE_INLINE
		explicit arg_forwarder(const Arg& a) : arg(a) { }
	};


	template<class Arg>
	LMAT_ENSURE_INLINE
	inline arg_forwarder<cref_arg_t, Arg> cref_arg(const Arg& arg)
	{
		return arg_forwarder<cref_arg_t, Arg>(arg);
	}

	template<class Arg>
	LMAT_ENSURE_INLINE
	inline arg_forwarder<ref_arg_t, Arg> ref_arg(Arg& arg)
	{
		return arg_forwarder<ref_arg_t, Arg>(arg);
	}


	template<class Arg>
	LMAT_ENSURE_INLINE
	inline arg_forwarder<cpy_arg_t, Arg> copy_arg(const Arg& arg)
	{
		return arg_forwarder<cpy_arg_t, Arg>(arg);
	}


	// Holders

	template<class Arg>
	class arg_holder<cref_arg_t, Arg>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit arg_holder(arg_forwarder<cref_arg_t, Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const Arg& get() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	template<class Arg>
	class arg_holder<ref_arg_t, Arg>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit arg_holder(arg_forwarder<ref_arg_t, Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		Arg& get() const
		{
			return m_arg;
		}

	private:
		const Arg& m_arg;
	};


	template<class Arg>
	class arg_holder<cpy_arg_t, Arg>
	{
	public:
		LMAT_ENSURE_INLINE
		explicit arg_holder(arg_forwarder<cpy_arg_t, Arg> fwd)
		: m_arg(fwd.arg)
		{
		}

		LMAT_ENSURE_INLINE
		const Arg& get() const
		{
			return m_arg;
		}

	private:
		Arg m_arg;
	};


	/********************************************
	 *
	 *  Argument list classes
	 *
	 ********************************************/

	// argument lists

	template<class QArgList> class qarg_list;

	template<class QArg1>
	class qarg_list< LMAT_TYPELIST_1(QArg1) >
	{
	public:
		LMAT_ENSURE_INLINE
		qarg_list(typename QArg1::forwarder arg1_fwd)
		: m_arg1_holder(arg1_fwd) { }

		LMAT_ENSURE_INLINE
		typename QArg1::reference arg1() const
		{
			return m_arg1_holder.get();
		}

	private:
		typename QArg1::holder m_arg1_holder;
	};


	template<class QArg1, class QArg2>
	class qarg_list< LMAT_TYPELIST_2(QArg1, QArg2) >
	{
	public:
		LMAT_ENSURE_INLINE
		qarg_list(
				typename QArg1::forwarder arg1_fwd,
				typename QArg2::forwarder arg2_fwd)
		: m_arg1_holder(arg1_fwd)
		, m_arg2_holder(arg2_fwd) { }

		LMAT_ENSURE_INLINE
		typename QArg1::reference arg1() const
		{
			return m_arg1_holder.get();
		}

		LMAT_ENSURE_INLINE
		typename QArg2::reference arg2() const
		{
			return m_arg2_holder.get();
		}

	private:
		typename QArg1::holder m_arg1_holder;
		typename QArg2::holder m_arg2_holder;
	};


	template<class QArg1, class QArg2, class QArg3>
	class qarg_list< LMAT_TYPELIST_3(QArg1, QArg2, QArg3) >
	{
	public:
		LMAT_ENSURE_INLINE
		qarg_list(
				typename QArg1::forwarder arg1_fwd,
				typename QArg2::forwarder arg2_fwd,
				typename QArg3::forwarder arg3_fwd)
		: m_arg1_holder(arg1_fwd)
		, m_arg2_holder(arg2_fwd)
		, m_arg3_holder(arg3_fwd) { }

		LMAT_ENSURE_INLINE
		typename QArg1::reference arg1() const
		{
			return m_arg1_holder.get();
		}

		LMAT_ENSURE_INLINE
		typename QArg2::reference arg2() const
		{
			return m_arg2_holder.get();
		}

		LMAT_ENSURE_INLINE
		typename QArg3::reference arg3() const
		{
			return m_arg3_holder.get();
		}

	private:
		typename QArg1::holder m_arg1_holder;
		typename QArg2::holder m_arg2_holder;
		typename QArg3::holder m_arg3_holder;
	};


	namespace meta
	{
		template<class Q>
		struct argument_of
		{
			typedef typename Q::argument type;
		};

		template<class Expr> struct num_args_;

		template<class QLst>
		struct num_args_<qarg_list<QLst> >
		{
			static const int value = length_<QLst>::value;
		};

		template<class Expr> struct get_argtype_list_;

		template<class QLst>
		struct get_argtype_list_<qarg_list<QLst> >
		{
			typedef typename transform_<argument_of, QLst>::type type;
		};
	}


	/********************************************
	 *
	 *  Expression construction
	 *
	 ********************************************/

	template<typename Spec, class ArgList> class expr_verifier;
	template<typename Spec, class QLst> class expr_map;

	template<typename Spec, typename Arg1_HP, class Arg1>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<
		expr_verifier<Spec, LMAT_TYPELIST_1(Arg1)>,
		typename expr_map<Spec,
			LMAT_TYPELIST_1( LMAT_QARG(Arg1_HP, Arg1) )>::type
	>::type
	make_expr(arg_forwarder<Arg1_HP, Arg1> arg1_fwd)
	{
		return expr_map<Spec,
				LMAT_TYPELIST_1( LMAT_QARG(Arg1_HP, Arg1) )
			>::get(arg1_fwd);
	}


	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<
		expr_verifier<Spec, LMAT_TYPELIST_2(Arg1, Arg2)>,
		typename expr_map<Spec,
			LMAT_TYPELIST_2(
					LMAT_QARG(Arg1_HP, Arg1),
					LMAT_QARG(Arg2_HP, Arg2)) >::type
	>::type
	make_expr(
			arg_forwarder<Arg1_HP, Arg1> arg1_fwd,
			arg_forwarder<Arg2_HP, Arg2> arg2_fwd)
	{
		return expr_map<Spec,
				LMAT_TYPELIST_2(
						LMAT_QARG(Arg1_HP, Arg1),
						LMAT_QARG(Arg2_HP, Arg2) )
			>::get(arg1_fwd, arg2_fwd);
	}


	template<typename Spec,
		typename Arg1_HP, class Arg1,
		typename Arg2_HP, class Arg2,
		typename Arg3_HP, class Arg3>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<
		expr_verifier<Spec, LMAT_TYPELIST_3(Arg1, Arg2, Arg3)>,
		typename expr_map<Spec,
			LMAT_TYPELIST_3(
					LMAT_QARG(Arg1_HP, Arg1),
					LMAT_QARG(Arg2_HP, Arg2),
					LMAT_QARG(Arg3_HP, Arg3)) >::type
	>::type
	make_expr(
			arg_forwarder<Arg1_HP, Arg1> arg1_fwd,
			arg_forwarder<Arg2_HP, Arg2> arg2_fwd,
			arg_forwarder<Arg3_HP, Arg3> arg3_fwd)
	{
		return expr_map<Spec,
				LMAT_TYPELIST_3(
						LMAT_QARG(Arg1_HP, Arg1),
						LMAT_QARG(Arg2_HP, Arg2),
						LMAT_QARG(Arg3_HP, Arg3))
			>::get(arg1_fwd, arg2_fwd, arg3_fwd);
	}

}

#endif /* EXPR_BASE_H_ */
