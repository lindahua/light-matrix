/**
 * @file matrix_print.h
 *
 * Functions to print matrices
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATRIX_PRINT_H_
#define LIGHTMAT_MATRIX_PRINT_H_

#include <light_mat/matrix/matrix_concepts.h>
#include <cstdio>


#define LMAT_DISP(a) ::lmat::disp_mat(#a, a)

#define LMAT_DEFINE_DEFAULT_MATPRINT_FMT( T, fmt ) \
	template<> struct default_matprint_fmt<T> { \
		static const char *get() { return fmt; } };


namespace lmat
{
	template<typename T, class Mat>
	inline void printf_mat(const char *fmt, const IRegularMatrix<Mat, T>& X,
			const char *pre_line=nullptr, const char *delim="\n")
	{
		index_t m = X.nrows();
		index_t n = X.ncolumns();

		for (index_t i = 0; i < m; ++i)
		{
			if (pre_line) std::printf("%s", pre_line);
			for (index_t j = 0; j < n; ++j)
			{
				std::printf(fmt, X.elem(i, j));
			}
			if (delim) std::printf("%s", delim);
		}
	}


	template<typename T> struct default_matprint_fmt;

	LMAT_DEFINE_DEFAULT_MATPRINT_FMT( double, "%10.4g" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT( float,  "%10.4g" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT(  int32_t,  "%8d" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT( uint32_t,  "%8u" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT(  int16_t,  "%6d" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT( uint16_t,  "%6u" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT(   int8_t,  "%4d" )
	LMAT_DEFINE_DEFAULT_MATPRINT_FMT(  uint8_t,  "%4u" )

	template<typename T, class Mat>
	const char* get_default_matprint_fmt(const IRegularMatrix<Mat, T>& X)
	{
		return default_matprint_fmt<T>::get();
	}

	template<typename T, class Mat>
	inline void printf_mat(const IRegularMatrix<Mat, T>& X)
	{
		printf_mat(get_default_matprint_fmt(X), X);
	}

	template<typename T, class Mat>
	inline void disp_mat(const char *name, const IRegularMatrix<Mat, T>& X)
	{
		std::printf("%s = \n", name);
		printf_mat(get_default_matprint_fmt(X), X, "    ");
		std::printf("\n");
	}
}

#endif /* MATRIX_PRINT_H_ */
