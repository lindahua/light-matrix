/**
 * @file linalg_fwd.h
 *
 * @brief Forward declarations for Linear algebra
 *
 * @author Dahua Lin
 */

#ifndef LINALG_FWD_H_
#define LINALG_FWD_H_

#include <light_mat/matrix/matrix_properties.h>

namespace lmat
{
	namespace blas
	{
		struct trs
		{
			const char uplo;
			const char trans;
			const char diag;

			LMAT_ENSURE_INLINE
			trs(char uplo_)
			: uplo(uplo_), trans('N'), diag('N') { }

			LMAT_ENSURE_INLINE
			trs(char uplo_, char trans_)
			: uplo(uplo_), trans(trans_), diag('N') { }

			LMAT_ENSURE_INLINE
			trs(char uplo_, char trans_, char diag_)
			: uplo(uplo_), trans(trans_), diag(diag_) { }
		};
	}


}


#endif /* LINALG_FWD_H_ */
