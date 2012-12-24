/**
 * @file linalg_fwd.h
 *
 * @brief Forward declarations for Linear algebra
 *
 * @author Dahua Lin
 */

#ifndef LINALG_FWD_H_
#define LINALG_FWD_H_

#include <light_mat/config/config.h>
#include <light_mat/matrix/matrix_classes.h>

#ifdef LMAT_BLAS_ILP64
typedef long long blas_int;
#else
typedef int blas_int;
#endif



#define LMAT_BLAS_NAME(name) name
#define LMAT_LAPACK_NAME(name) LMAT_BLAS_NAME(name)

#define LMAT_CHECK_WHOLE_CONT(Ty) static_assert( meta::is_continuous<Ty>::value, #Ty " must be continuous.");
#define LMAT_CHECK_PERCOL_CONT(Ty) static_assert( meta::is_percol_continuous<Ty>::value, #Ty " must be percol continuous.");

namespace lmat
{

	namespace meta
	{
		template<class Mat>
		struct sq_dim
		{
			static const int value = meta::common_dim<
					meta::nrows<Mat>::value, meta::ncols<Mat>::value>::value;
		};
	}


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
