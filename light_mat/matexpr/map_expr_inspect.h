/**
 * @file map_expr_inspect.h
 *
 * @brief Inspection of map expression
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAP_EXPR_INSPECT_H_
#define LIGHTMAT_MAP_EXPR_INSPECT_H_

#include <light_mat/mateval/map_expr.h>
#include <light_mat/matrix/matrix_inspect.h>

//#define LMAT_DUMP_FUN

#define _LMAT_DEFINE_FTAG_NAME(T) \
	template<> struct type_name<ftags::T> { \
		static std::string get() { return #T; } };

namespace lmat
{

	/********************************************
	 *
	 *  define names for fun tags
	 *
	 ********************************************/

	// arithmetic

	_LMAT_DEFINE_FTAG_NAME( add_ )
	_LMAT_DEFINE_FTAG_NAME( sub_ )
	_LMAT_DEFINE_FTAG_NAME( mul_ )
	_LMAT_DEFINE_FTAG_NAME( div_ )
	_LMAT_DEFINE_FTAG_NAME( neg_ )

	_LMAT_DEFINE_FTAG_NAME( abs_ )
	_LMAT_DEFINE_FTAG_NAME( sqr_ )
	_LMAT_DEFINE_FTAG_NAME( cube_ )

	_LMAT_DEFINE_FTAG_NAME( max_ )
	_LMAT_DEFINE_FTAG_NAME( min_ )
	_LMAT_DEFINE_FTAG_NAME( clamp_ )

	// comparison

	_LMAT_DEFINE_FTAG_NAME( eq_ )
	_LMAT_DEFINE_FTAG_NAME( ne_ )
	_LMAT_DEFINE_FTAG_NAME( ge_ )
	_LMAT_DEFINE_FTAG_NAME( gt_ )
	_LMAT_DEFINE_FTAG_NAME( le_ )
	_LMAT_DEFINE_FTAG_NAME( lt_ )

	// logical

	_LMAT_DEFINE_FTAG_NAME( logical_not_ )
	_LMAT_DEFINE_FTAG_NAME( logical_and_ )
	_LMAT_DEFINE_FTAG_NAME( logical_or_ )
	_LMAT_DEFINE_FTAG_NAME( logical_eq_ )
	_LMAT_DEFINE_FTAG_NAME( logical_ne_ )


	// real math

	_LMAT_DEFINE_FTAG_NAME( fma_ )
	_LMAT_DEFINE_FTAG_NAME( rcp_ )
	_LMAT_DEFINE_FTAG_NAME( sqrt_ )
	_LMAT_DEFINE_FTAG_NAME( rsqrt_ )
	_LMAT_DEFINE_FTAG_NAME( pow_ )

	_LMAT_DEFINE_FTAG_NAME( floor_ )
	_LMAT_DEFINE_FTAG_NAME( ceil_ )

	_LMAT_DEFINE_FTAG_NAME( exp_ )
	_LMAT_DEFINE_FTAG_NAME( log_ )
	_LMAT_DEFINE_FTAG_NAME( log10_ )
	_LMAT_DEFINE_FTAG_NAME( xlogx_ )
	_LMAT_DEFINE_FTAG_NAME( xlogy_ )

	_LMAT_DEFINE_FTAG_NAME( sin_ )
	_LMAT_DEFINE_FTAG_NAME( cos_ )
	_LMAT_DEFINE_FTAG_NAME( tan_ )

	_LMAT_DEFINE_FTAG_NAME( asin_ )
	_LMAT_DEFINE_FTAG_NAME( acos_ )
	_LMAT_DEFINE_FTAG_NAME( atan_ )
	_LMAT_DEFINE_FTAG_NAME( atan2_ )

	_LMAT_DEFINE_FTAG_NAME( sinh_ )
	_LMAT_DEFINE_FTAG_NAME( cosh_ )
	_LMAT_DEFINE_FTAG_NAME( tanh_ )

	// C++11 real math

	_LMAT_DEFINE_FTAG_NAME( cbrt_ )
	_LMAT_DEFINE_FTAG_NAME( hypot_ )

	_LMAT_DEFINE_FTAG_NAME( round_ )
	_LMAT_DEFINE_FTAG_NAME( trunc_ )

	_LMAT_DEFINE_FTAG_NAME( exp2_ )
	_LMAT_DEFINE_FTAG_NAME( log2_ )
	_LMAT_DEFINE_FTAG_NAME( expm1_ )
	_LMAT_DEFINE_FTAG_NAME( log1p_ )

	_LMAT_DEFINE_FTAG_NAME( asinh_ )
	_LMAT_DEFINE_FTAG_NAME( acosh_ )
	_LMAT_DEFINE_FTAG_NAME( atanh_ )

	_LMAT_DEFINE_FTAG_NAME( erf_ )
	_LMAT_DEFINE_FTAG_NAME( erfc_ )
	_LMAT_DEFINE_FTAG_NAME( lgamma_ )
	_LMAT_DEFINE_FTAG_NAME( tgamma_ )

	// numeric predicates

	_LMAT_DEFINE_FTAG_NAME( signbit_ )
	_LMAT_DEFINE_FTAG_NAME( isfinite_ )
	_LMAT_DEFINE_FTAG_NAME( isinf_ )
	_LMAT_DEFINE_FTAG_NAME( isnan_ )

	_LMAT_DEFINE_FTAG_NAME( cond_ )


	/********************************************
	 *
	 *  general devices
	 *
	 ********************************************/

	template<typename FTag, typename... Args>
	struct expr_name<map_expr<FTag, Args...> >
	{
		static std::string get()
		{
			return std::string("map(") + type_name<FTag>::get() + ")";
		}
	};

	template<typename FTag, typename... Args>
	inline void dump_top_(std::ostream& out, const map_expr<FTag, Args...>& expr, int indent)
	{
		dump_top(out, expr, indent);

#ifdef LMAT_DUMP_FUN
		typedef typename fun_map<FTag, typename internal::arg_value_type<Args>::type...>::type fun_t;

		dump_indent(out, indent+1);
		out << "Fun: " << type_name<fun_t>::get()
			<< " is_simdizable = " << is_simdizable<fun_t, default_simd_kind>::value << "\n";
#endif
	}


	template<typename FTag, typename A1>
	inline void dump_expr(std::ostream& out, const map_expr<FTag, A1>& expr, int indent=0)
	{
		dump_top_(out, expr, indent);
		dump_expr(out, expr.arg1(), indent+1);
	}

	template<typename FTag, typename A1, typename A2>
	inline void dump_expr(std::ostream& out, const map_expr<FTag, A1, A2>& expr, int indent=0)
	{
		dump_top_(out, expr, indent);
		dump_expr(out, expr.arg1(), indent+1);
		dump_expr(out, expr.arg2(), indent+1);
	}

	template<typename FTag, typename A1, typename A2, typename A3>
	inline void dump_expr(std::ostream& out, const map_expr<FTag, A1, A2, A3>& expr, int indent=0)
	{
		dump_top_(out, expr, indent);
		dump_expr(out, expr.arg1(), indent+1);
		dump_expr(out, expr.arg2(), indent+1);
		dump_expr(out, expr.arg3(), indent+1);
	}

}

#endif
