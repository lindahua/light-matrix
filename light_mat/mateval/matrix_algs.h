/**
 * @file matrix_algs.h
 *
 * @brief Generic algorithms on matrices
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MATRIX_ALGS_H_
#define LIGHTMAT_MATRIX_ALGS_H_

#include <light_mat/matrix/matrix_classes.h>
#include "internal/matrix_algs_internal.h"

#include <functional>
#include <algorithm>


namespace lmat
{

	/********************************************
	 *
	 *  counting
	 *
	 ********************************************/

	template<class A, typename T>
	inline size_t count(const IEWiseMatrix<A, T>& a)
	{
		return internal::count_impl<supports_linear_macc<A>::value>::run(a.derived());
	}

	template<class A, typename T, class D, typename TD>
	inline void colwise_count(const IEWiseMatrix<A, T>& a, IRegularMatrix<D, TD>& dmat)
	{
		D& d = dmat.derived();
		const index_t m = a.nrows();
		const index_t n = a.ncolumns();

		LMAT_CHECK_DIMS( dmat.nelems() == n )
		auto rd = make_multicol_accessor(atags::scalar(), in_(a.derived()));

		for (index_t j = 0; j < n; ++j)
		{
			dmat[j] = internal::count_by_reader(m, rd.col(j));
		}
	}


}

#endif /* MATRIX_ALGS_H_ */
