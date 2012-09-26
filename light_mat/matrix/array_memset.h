/*
 * @file array_memset.h
 *
 * @brief Array memory setting
 *
 * @author Dahua Lin
 */

#ifndef ARRAY_MEMSET_H_
#define ARRAY_MEMSET_H_

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/common/memory.h>

namespace lmat
{
	template<typename T, class Arr, class Setter>
	inline void set_array_memory(IDenseMatrix<Arr, T>& arr, IMemorySetter<Setter, T> setter)
	{
		if (has_continuous_layout(arr))
		{
			setter.set(arr.nelems(), arr.ptr_data());
		}
		else
		{
			T *p = arr.ptr_data();
			const index_t m = arr.nrows();
			const index_t n = arr.ncolumns();
			const index_t ldim = arr.lead_dim();

			if (m == 1)
			{
				for (index_t j = 0; j < n; ++j, p += ldim)
					setter.set(*p);
			}
			else
			{
				for (index_t j = 0; j < n; ++j, p += ldim)
					setter.set(m, p);
			}
		}
	}

}

#endif /* ARRAY_MEMSET_H_ */
