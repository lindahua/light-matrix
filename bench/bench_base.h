/**
 * @file bench_base.h
 *
 * @brief The basis for benchmark
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_BENCH_BASE_H_
#define LIGHTMAT_BENCH_BASE_H_

#include <light_mat/matrix/matrix_classes.h>
#include <light_test/benchmark.h>
#include <light_test/std_bench_mon.h>

namespace lmat { namespace bench {

	LMAT_ENSURE_INLINE
	inline double rand_unif()
	{
		return double(std::rand()) / RAND_MAX;
	}

	template<typename T, class Mat>
	inline void fill_rand(IRegularMatrix<Mat, T>& mat)
	{
		const index_t m = mat.nrows();
		const index_t n = mat.ncolumns();

		for (index_t j = 0; j < n; ++j)
		{
			for (index_t i = 0; i < m; ++i)
			{
				mat(i, j) = (T)rand_unif();
			}
		}
	}

} }

#endif /* BENCH_BASE_H_ */
