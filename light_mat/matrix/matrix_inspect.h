/**
 * @file matrix_inspect.h
 *
 * @brief Inspection of matrix expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_INSPECT_H_
#define LIGHTMAT_MATRIX_INSPECT_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/mateval/macc_policy.h>
#include <string>
#include <iostream>
#include <typeinfo>


#define _LMAT_DEFINE_TYPE_NAME(T) \
	template<> struct type_name<T> { \
		static std::string get() { return #T; } };

#define _LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( TC ) \
	template<typename T, int M, int N> \
	struct expr_name< TC<T, M, N> > { \
		static std::string get() { return #TC; } };

namespace lmat
{
	// type names

	template<typename T>
	struct type_name
	{
		static std::string get()
		{
			return typeid(T).name();
		}
	};

	_LMAT_DEFINE_TYPE_NAME(float)
	_LMAT_DEFINE_TYPE_NAME(double)

	_LMAT_DEFINE_TYPE_NAME(int8_t)
	_LMAT_DEFINE_TYPE_NAME(uint8_t)
	_LMAT_DEFINE_TYPE_NAME(int16_t)
	_LMAT_DEFINE_TYPE_NAME(uint16_t)
	_LMAT_DEFINE_TYPE_NAME(int32_t)
	_LMAT_DEFINE_TYPE_NAME(uint32_t)
	_LMAT_DEFINE_TYPE_NAME(int64_t)
	_LMAT_DEFINE_TYPE_NAME(uint64_t)

	template<typename T>
	struct type_name<mask_t<T> >
	{
		static std::string get()
		{
			return std::string("mask_t<") + type_name<T>::get() + ">";
		}
	};


	// expression names

	template<class Expr>
	struct expr_name
	{
		static std::string get()
		{
			return "unknown_expr";
		}
	};

	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( dense_matrix )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( cref_matrix )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( ref_matrix )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( cref_block )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( ref_block )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( cref_grid )
	_LMAT_DEFINE_EXPR_NAME_FOR_REGMAT( ref_grid )

	// dump expression

	template<class Expr, typename T>
	inline std::string get_expr_name(const IEWiseMatrix<Expr, T>& )
	{
		return expr_name<Expr>::get();
	}

	inline void dump_indent(std::ostream& out, int indent)
	{
		for (int i = 0; i < indent; ++i) out << "    ";
	}

	template<class Expr, typename T>
	inline void dump_top(std::ostream& out, const IEWiseMatrix<Expr, T>& expr, int indent=0)
	{
		dump_indent(out, indent);

		typedef preferred_macc_policy<Expr> pmap;
		const char *acc_pat = pmap::prefer_linear ? "linear" : "percol";
		const char* acc_u = pmap::prefer_simd ? "simd" : "scalar";

		out << get_expr_name(expr)
			<< "["
				<< type_name<T>::get() << ' '
				<< meta::nrows<Expr>::value << " x " << meta::ncols<Expr>::value
			<< "]: "
			<< "m = " << expr.nrows() << ", "
			<< "n = " << expr.ncolumns() << ", "
			<< "macc = (" << acc_pat << ", " << acc_u << ")\n";
	}

	template<class Expr, typename T>
	inline void dump_expr(std::ostream& out, const IEWiseMatrix<Expr, T>& expr, int indent=0)
	{
		dump_top(out, expr, indent);
	}


}

#endif
