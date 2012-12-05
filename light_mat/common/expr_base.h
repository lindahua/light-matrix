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
	 *  Tied forwarder
	 *
	 ********************************************/

	template<class QLst> struct tied_forwarder;

	template<class QArg1>
	struct tied_forwarder< LMAT_TYPELIST_1(QArg1) >
	{
		tied_forwarder(typename QArg1::reference a1)
		: arg1_fwd(a1) { }

		typename QArg1::forwarder arg1_fwd;
	};

	template<class QArg1, class QArg2>
	struct tied_forwarder< LMAT_TYPELIST_2(QArg1, QArg2) >
	{
		tied_forwarder(typename QArg1::reference a1, typename QArg2::reference a2)
		: arg1_fwd(a1), arg2_fwd(a2) { }

		typename QArg1::forwarder arg1_fwd;
		typename QArg2::forwarder arg2_fwd;
	};

	template<class QArg1, class QArg2, class QArg3>
	struct tied_forwarder< LMAT_TYPELIST_3(QArg1, QArg2, QArg3) >
	{
		tied_forwarder(
				typename QArg1::reference a1,
				typename QArg2::reference a2,
				typename QArg3::reference a3)
		: arg1_fwd(a1), arg2_fwd(a2), arg3_fwd(a3) { }

		typename QArg1::forwarder arg1_fwd;
		typename QArg2::forwarder arg2_fwd;
		typename QArg3::forwarder arg3_fwd;
	};

	template<class QArg1, class QArg2, class QArg3, class QArg4>
	struct tied_forwarder< LMAT_TYPELIST_4(QArg1, QArg2, QArg3, QArg4) >
	{
		tied_forwarder(
				typename QArg1::reference a1,
				typename QArg2::reference a2,
				typename QArg3::reference a3,
				typename QArg4::reference a4 )
		: arg1_fwd(a1), arg2_fwd(a2), arg3_fwd(a3), arg4_fwd(a4) { }

		typename QArg1::forwarder arg1_fwd;
		typename QArg2::forwarder arg2_fwd;
		typename QArg3::forwarder arg3_fwd;
		typename QArg4::forwarder arg4_fwd;
	};


	/********************************************
	 *
	 *  Argument list classes
	 *
	 ********************************************/

	template<class QLst> class qarg_list;

	template<class QArg1>
	class qarg_list< LMAT_TYPELIST_1(QArg1) >
	{
	public:
		typedef tied_forwarder< LMAT_TYPELIST_1(QArg1) > tied_forwarder_t;

		LMAT_ENSURE_INLINE
		qarg_list(tied_forwarder_t fwd)
		: m_arg1_holder(fwd.arg1_fwd) { }

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
		typedef tied_forwarder< LMAT_TYPELIST_2(QArg1, QArg2) > tied_forwarder_t;

		LMAT_ENSURE_INLINE
		qarg_list( tied_forwarder_t fwd )
		: m_arg1_holder(fwd.arg1_fwd)
		, m_arg2_holder(fwd.arg2_fwd) { }

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
		typedef tied_forwarder< LMAT_TYPELIST_3(QArg1, QArg2, QArg3) > tied_forwarder_t;

		LMAT_ENSURE_INLINE
		qarg_list( tied_forwarder_t fwd )
		: m_arg1_holder(fwd.arg1_fwd)
		, m_arg2_holder(fwd.arg2_fwd)
		, m_arg3_holder(fwd.arg3_fwd) { }

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


	template<class QArg1, class QArg2, class QArg3, class QArg4>
	class qarg_list< LMAT_TYPELIST_4(QArg1, QArg2, QArg3, QArg4) >
	{
	public:
		typedef tied_forwarder< LMAT_TYPELIST_4(QArg1, QArg2, QArg3, QArg4) > tied_forwarder_t;

		LMAT_ENSURE_INLINE
		qarg_list( tied_forwarder_t fwd )
		: m_arg1_holder(fwd.arg1_fwd)
		, m_arg2_holder(fwd.arg2_fwd)
		, m_arg3_holder(fwd.arg3_fwd)
		, m_arg4_holder(fwd.arg4_fwd) { }

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

		LMAT_ENSURE_INLINE
		typename QArg4::reference arg4() const
		{
			return m_arg4_holder.get();
		}

	private:
		typename QArg1::holder m_arg1_holder;
		typename QArg2::holder m_arg2_holder;
		typename QArg3::holder m_arg3_holder;
		typename QArg4::holder m_arg4_holder;
	};


	namespace meta
	{
		template<class Q>
		struct argument_of
		{
			typedef typename Q::argument type;
		};

		template<class QLst>
		struct num_args_
		{
			static const int value = length_<QLst>::value;
		};

		template<class QLst>
		struct argtype_list_
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

	template<typename Spec, class QLst>
	LMAT_ENSURE_INLINE
	typename meta::enable_if<
		expr_verifier<Spec, typename meta::argtype_list_<QLst>::type >,
		typename expr_map<Spec, QLst>::type
	>::type
	make_expr(const Spec& spec, const tied_forwarder<QLst>& tfwd)
	{
		return expr_map<Spec, QLst>::get(spec, tfwd);
	}

}

#endif /* EXPR_BASE_H_ */
