/**
 * @file rand_accessors.h
 *
 * @brief Matrix accessors for random expressions
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_RAND_ACCESSORS_H_
#define LIGHTMAT_RAND_ACCESSORS_H_

#include <light_mat/mateval/multicol_accessors.h>
#include <light_mat/random/distr_fwd.h>

namespace lmat
{
	// forward declarations

	template<typename RStream, typename Distr, typename U> class rand_vec_reader;
	template<typename RStream, typename Distr, typename U> class rand_multicol_reader;

	/********************************************
	 *
	 *  Vector readers
	 *
	 ********************************************/

	template<typename RStream, typename Distr>
	class rand_vec_reader<RStream, Distr, scalar_> : public scalar_vec_accessor_base
	{
		typedef scalar_ atag;
		typedef typename Distr::result_type result_t;

	public:
		LMAT_ENSURE_INLINE
		rand_vec_reader(RStream& rstream, const Distr& distr)
		: m_rstream(rstream), m_distr(distr)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_distr(m_rstream);
		}

	private:
		RStream& m_rstream;
		const Distr& m_distr;
	};

	template<typename RStream, typename Distr, typename Kind>
	class rand_vec_reader<RStream, Distr, simd_<Kind>> : public simd_vec_accessor_base
	{
		typedef simd_<Kind> atag;

		typedef typename Distr::result_type result_t;
		typedef typename simdize_map<Distr, Kind>::type simd_distr_t;
		typedef typename simd_distr_t::result_type pack_t;

	public:
		LMAT_ENSURE_INLINE
		rand_vec_reader(RStream& rstream, const Distr& distr)
		: m_rstream(rstream), m_distr(distr)
		, m_distr_simd(simdize_map<Distr, Kind>::get(distr))
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_distr(m_rstream);
		}

		LMAT_ENSURE_INLINE
		pack_t pack(index_t i) const
		{
			return m_distr_simd(m_rstream);
		}

	private:
		RStream& m_rstream;
		const Distr& m_distr;
		simd_distr_t m_distr_simd;
	};


	/********************************************
	 *
	 *  Multi-column readers
	 *
	 ********************************************/

	template<typename RStream, typename Distr, typename U>
	class rand_multicol_reader : public multicol_accessor_base
	{
	public:
		typedef rand_vec_reader<RStream, Distr, U> col_accessor_type;

		LMAT_ENSURE_INLINE
		rand_multicol_reader(RStream& rstream, const Distr& distr)
		: m_rstream(rstream), m_distr(distr) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t) const
		{
			return col_accessor_type(m_rstream, m_distr);
		}

	private:
		RStream& m_rstream;
		const Distr& m_distr;
	};
}

#endif
