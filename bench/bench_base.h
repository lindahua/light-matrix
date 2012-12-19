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
#include <light_mat/common/timer.h>
#include <cstdio>
#include <cstdlib>

namespace lmat { namespace bench {


	enum performance_unit
	{
		PUNIT_NONE = 0,
		PUNIT_KPS = 1,
		PUNIT_MPS = 2,
		PUNIT_GPS = 3
	};


	inline double et_to_speed(performance_unit u, uint64_t N, double et)
	{
		double c = 0;

		switch (u)
		{
		case PUNIT_KPS:
			c = 1.0e-3;
			break;

		case PUNIT_MPS:
			c = 1.0e-6;
			break;

		case PUNIT_GPS:
			c = 1.0e-9;
			break;

		case PUNIT_NONE:
			break;
		}

		return double(N) * c / et;
	}


	inline const char *punit_to_name(performance_unit u)
	{
		switch (u)
		{
		case PUNIT_KPS:
			return "KPS";

		case PUNIT_MPS:
			return "MPS";

		case PUNIT_GPS:
			return "GPS";

		case PUNIT_NONE:
			break;
		}

		return "";
	}


	class std_bench_mon
	{
	public:
		std_bench_mon(performance_unit perfu)
		: m_perfu(perfu) { }

		template<class Job>
		void report(const Job& job, unsigned int nrun, double et)
		{
			std::printf("%20s:  %8.4f sec / %9u run", job.name(), et, nrun);

			uint64_t N = uint64_t(job.size()) * uint64_t(nrun);

			if (m_perfu != PUNIT_NONE)
			{
				std::printf("  |  %9.3f %s\n", et_to_speed(m_perfu, N, et), punit_to_name(m_perfu));
			}
			else
			{
				std::printf("\n");
			}
		}

	private:
		performance_unit m_perfu;
	};

	template<class Job, class Mon>
	inline void time_job(Job& job, unsigned int nwarm, unsigned int nrun, Mon& mon)
	{
		for (unsigned int i = 0; i < nwarm; ++i) job.run();

		timer tm(true);
		for (unsigned int i = 0; i < nrun; ++i) job.run();
		double et = tm.elapsed_secs();

		mon.report(job, nrun, et);
	}

	template<class Job>
	inline void time_job(Job& job, unsigned int nwarm, unsigned int nrun, performance_unit perfu)
	{
		std_bench_mon mon(perfu);
		time_job(job, nwarm, nrun, mon);
	}


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
