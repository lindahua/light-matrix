/**
 * @file timer.h
 *
 * @brief High quality timer for benchmark
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_TIMER_H_
#define LIGHTMAT_TIMER_H_

#include <light_mat/common/prim_types.h>

#include <time.h>
#include <sys/time.h>

#ifdef __MACH__
#include <mach/mach_time.h>
#endif

namespace lmat
{

#ifdef __MACH__

	class timer
	{
	public:
		LMAT_ENSURE_INLINE
		explicit timer( bool to_start = false )
		{
			::mach_timebase_info(&m_baseinfo);
			if (to_start) start();
		}

		LMAT_ENSURE_INLINE
		void start()
		{
			m_start_t = ::mach_absolute_time();
		}

		LMAT_ENSURE_INLINE
		double elapsed_secs() const
		{
			uint64_t t = ::mach_absolute_time();
			return 1.0e-9 * double(to_nanosecs(t));
		}

		LMAT_ENSURE_INLINE
		double elapsed_msecs() const
		{
			uint64_t t = ::mach_absolute_time();
			return 1.0e-6 * double(to_nanosecs(t));
		}

		LMAT_ENSURE_INLINE
		double elapsed_usecs() const
		{
			uint64_t t = ::mach_absolute_time();
			return 1.0e-3 * double(to_nanosecs(t));
		}

	private:
		LMAT_ENSURE_INLINE
		uint64_t to_nanosecs(const uint64_t& t) const
		{
			return (m_baseinfo.numer * (t - m_start_t)) / m_baseinfo.denom;
		}

	private:
		uint64_t m_start_t;
		mach_timebase_info_data_t m_baseinfo;
	};

#else

	class timer
	{
	public:
		LMAT_ENSURE_INLINE
		explicit timer( bool to_start = false )
		{
			if (to_start) start();
		}

		LMAT_ENSURE_INLINE
		void start()
		{
			::clock_gettime(CLOCK_REALTIME, &_start_tp);
		}

		LMAT_ENSURE_INLINE
		double elapsed_secs() const
		{
			timespec_t t;
			::clock_gettime(CLOCK_REALTIME, &t);
			return 1.0e-9 * to_nanosecs(t);
		}

		LMAT_ENSURE_INLINE
		double elapsed_msecs() const
		{
			timespec_t t;
			::clock_gettime(CLOCK_REALTIME, &t);
			return 1.0e-6 * to_nanosecs(t);
		}

		LMAT_ENSURE_INLINE
		double elapsed_usecs() const
		{
			timespec_t t;
			::clock_gettime(CLOCK_REALTIME, &t);
			return 1.0e-3 * to_nanosecs(t);
		}

	private:
		LMAT_ENSURE_INLINE
		double to_nanosecs(const timespec_t& t) const
		{
			return double(tp.tv_sec - _start_tp.tv_sec) * 1.0e9 +
					double(tp.tv_nsec - _start_tp.tv_nsec);
		}

	private:
		timespec_t _start_tp;
	};


#endif

}

#endif /* TIMER_H_ */
