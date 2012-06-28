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

namespace lmat
{
	template<typename T, class Mat>
	inline void printf_mat(const char *fmt, const IMatrixView<Mat, T>& X,
			const char *pre_line=LMAT_NULL, const char *delim="\n")
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
}

#endif /* MATRIX_PRINT_H_ */
