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

typedef blas_int lapack_int;

#define LMAT_BLAS_NAME(name) name
#define LMAT_LAPACK_NAME(name) LMAT_BLAS_NAME(name)

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


	namespace lapack
	{
		class lapack_failure : public std::exception
		{
		public:
			lapack_failure(const char *routine_name, int code)
			: m_routine(routine_name)
			, m_msg(std::string(routine_name) + " failed.")
			, m_errcode(code) { }

			int error_code() const throw()
			{
				return m_errcode;
			}

			virtual const char *what() const throw()
			{
				return m_msg.c_str();
			}

		private:
			std::string m_routine;
			std::string m_msg;
			int m_errcode;
		};
	}

}


#define LMAT_CHECK_LAPACK_INFO( fun, v ) if (v != 0) throw lmat::lapack::lapack_failure(#fun, (int)(v))


#endif /* LINALG_FWD_H_ */
