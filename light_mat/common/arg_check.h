/**
 * @file arg_check.h
 *
 * Argument checking
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_ARG_CHECK_H_
#define LIGHTMAT_ARG_CHECK_H_

#include <light_mat/common/lang_base.h>
#include <exception>

namespace lmat
{
	class invalid_operation : public std::exception
	{
	public:
		LMAT_ENSURE_INLINE
		invalid_operation(const char *msg)
		: m_msg(msg)
		{
		}

		virtual const char* what() const throw()
		{
			return m_msg;
		}

	private:
		const char *m_msg;
	};


	class invalid_argument : public std::exception
	{
	public:
		LMAT_ENSURE_INLINE
		invalid_argument(const char *msg)
		: m_msg(msg)
		{
		}

		virtual const char* what() const throw()
		{
			return m_msg;
		}

	private:
		const char *m_msg;
	};


	class out_of_range : public std::exception
	{
	public:
		out_of_range(const char *msg)
		: m_msg(msg)
		{
		}

		virtual const char* what() const throw()
		{
			return m_msg;
		}

	private:
		const char *m_msg;
	};


	// generic argument checking


	LMAT_ENSURE_INLINE inline void check_arg(bool cond, const char* message)
	{
		if (!cond)
		{
			throw invalid_argument(message);
		}
	}

	LMAT_ENSURE_INLINE inline void check_range(bool cond, const char *message)
	{
		if (!cond)
		{
			throw out_of_range(message);
		}
	}

	template<typename T>
	LMAT_ENSURE_INLINE inline const T& check_forward(const T& val, bool cond, const char *message)
	{
		check_arg(cond, message);
		return val;
	}


	template<bool Cond, typename T> struct conditional_enabled_value;

	template<typename T> struct conditional_enabled_value<true, T>
	{
		LMAT_ENSURE_INLINE
		static const T& get(const T& v, const char *msg)
		{
			return v;
		}
	};

	template<typename T> struct conditional_enabled_value<false, T>
	{
		LMAT_ENSURE_INLINE
		static const T& get(const T& v, const char *msg)
		{
			throw invalid_operation(msg);
		}
	};

}



#endif
