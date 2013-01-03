/**
 * @file lapack_fwd.h
 *
 * @brief Forward header for LAPACK functions/classes
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_LAPACK_FWD_H_
#define LIGHTMAT_LAPACK_FWD_H_

#include <light_mat/linalg/linalg_fwd.h>
#include "internal/linalg_aux.h"

#define LMAT_CALL_LAPACK(fun, params) \
	LMAT_LAPACK_NAME(fun) params; \
	if (info != 0) throw lmat::lapack::lapack_failure(#fun, (int)(info))

typedef blas_int lapack_int;

namespace lmat { namespace lapack {

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

} }

#endif /* LAPACK_FWD_H_ */
