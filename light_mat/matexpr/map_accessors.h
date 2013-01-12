/**
 * @file map_accessors.h
 *
 * @brief Function-mapping accessors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAP_ACCESSORS_H_
#define LIGHTMAT_MAP_ACCESSORS_H_

#include <light_mat/mateval/multicol_accessors.h>
#include <light_mat/math/functor_base.h>

namespace lmat
{
	// forward declarations

	template<typename Fun, typename U, typename... ArgReaders> class map_vec_reader;
	template<typename Fun, typename U, typename... ArgReaders> class map_multicol_reader;


	/********************************************
	 *
	 *  Vector readers
	 *
	 ********************************************/

	template<typename Fun, typename Rd1>
	class map_vec_reader<Fun, atags::scalar, Rd1> : public scalar_vec_accessor_base
	{
		typedef atags::scalar atag;
		typedef typename Fun::result_type result_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u, const Rd1& rd1)
		: m_fun(fun), m_rd1(rd1)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
	};

	template<typename Fun, typename Kind, typename Rd1>
	class map_vec_reader<Fun, atags::simd<Kind>, Rd1> : public simd_vec_accessor_base
	{
		typedef atags::simd<Kind> atag;

		typedef typename Fun::result_type result_t;
		typedef typename simdize_map<Fun, Kind>::type simd_fun_t;
		typedef typename simd_fun_t::result_type pack_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u, const Rd1& rd1)
		: m_fun(fun)
		, m_pkfun(simdize_map<Fun, Kind>::get(fun))
		, m_rd1(rd1)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i));
		}

		LMAT_ENSURE_INLINE
		pack_t pack(index_t i) const
		{
			return m_pkfun(m_rd1.pack(i));
		}

	private:
		Fun m_fun;
		simd_fun_t m_pkfun;
		Rd1 m_rd1;
	};


	template<typename Fun, typename Rd1, typename Rd2>
	class map_vec_reader<Fun, atags::scalar, Rd1, Rd2> : public scalar_vec_accessor_base
	{
		typedef atags::scalar atag;
		typedef typename Fun::result_type result_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u, const Rd1& rd1, const Rd2& rd2)
		: m_fun(fun), m_rd1(rd1), m_rd2(rd2)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i), m_rd2.scalar(i));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
		Rd2 m_rd2;
	};

	template<typename Fun, typename Kind, typename Rd1, typename Rd2>
	class map_vec_reader<Fun, atags::simd<Kind>, Rd1, Rd2> : public simd_vec_accessor_base
	{
		typedef atags::simd<Kind> atag;

		typedef typename Fun::result_type result_t;
		typedef typename simdize_map<Fun, Kind>::type simd_fun_t;
		typedef typename simd_fun_t::result_type pack_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u, const Rd1& rd1, const Rd2& rd2)
		: m_fun(fun)
		, m_pkfun(simdize_map<Fun, Kind>::get(fun))
		, m_rd1(rd1), m_rd2(rd2)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i), m_rd2.scalar(i));
		}

		LMAT_ENSURE_INLINE
		pack_t pack(index_t i) const
		{
			return m_pkfun(m_rd1.pack(i), m_rd2.pack(i));
		}

	private:
		Fun m_fun;
		simd_fun_t m_pkfun;
		Rd1 m_rd1;
		Rd2 m_rd2;
	};


	template<typename Fun, typename Rd1, typename Rd2, typename Rd3>
	class map_vec_reader<Fun, atags::scalar, Rd1, Rd2, Rd3> : public scalar_vec_accessor_base
	{
		typedef atags::scalar atag;
		typedef typename Fun::result_type result_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u,
				const Rd1& rd1, const Rd2& rd2, const Rd3& rd3)
		: m_fun(fun), m_rd1(rd1), m_rd2(rd2), m_rd3(rd3)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i), m_rd2.scalar(i), m_rd3.scalar(i));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
		Rd2 m_rd2;
		Rd3 m_rd3;
	};

	template<typename Fun, typename Kind, typename Rd1, typename Rd2, typename Rd3>
	class map_vec_reader<Fun, atags::simd<Kind>, Rd1, Rd2, Rd3> : public simd_vec_accessor_base
	{
		typedef atags::simd<Kind> atag;

		typedef typename Fun::result_type result_t;
		typedef typename simdize_map<Fun, Kind>::type simd_fun_t;
		typedef typename simd_fun_t::result_type pack_t;

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(const Fun& fun, atag u,
				const Rd1& rd1, const Rd2& rd2, const Rd3& rd3)
		: m_fun(fun)
		, m_pkfun(simdize_map<Fun, Kind>::get(fun))
		, m_rd1(rd1), m_rd2(rd2), m_rd3(rd3)
		{ }

		LMAT_ENSURE_INLINE
		result_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i), m_rd2.scalar(i), m_rd3.scalar(i));
		}

		LMAT_ENSURE_INLINE
		pack_t pack(index_t i) const
		{
			return m_pkfun(m_rd1.pack(i), m_rd2.pack(i), m_rd3.pack(i));
		}

	private:
		Fun m_fun;
		simd_fun_t m_pkfun;
		Rd1 m_rd1;
		Rd2 m_rd2;
		Rd3 m_rd3;
	};


	/********************************************
	 *
	 *  Multi-column readers
	 *
	 ********************************************/

	template<typename Fun, typename U, typename Rd1>
	class map_multicol_reader<Fun, U, Rd1> : public multicol_accessor_base
	{
		typedef typename Rd1::col_accessor_type col1_t;

	public:
		typedef map_vec_reader<Fun, U, col1_t> col_accessor_type;

		LMAT_ENSURE_INLINE
		map_multicol_reader(const Fun& fun, U, const Rd1& rd1)
		:m_fun(fun), m_rd1(rd1) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_fun, U(), m_rd1.col(j));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
	};


	template<typename Fun, typename U, typename Rd1, typename Rd2>
	class map_multicol_reader<Fun, U, Rd1, Rd2> : public multicol_accessor_base
	{
		typedef typename Rd1::col_accessor_type col1_t;
		typedef typename Rd2::col_accessor_type col2_t;

	public:
		typedef map_vec_reader<Fun, U, col1_t, col2_t> col_accessor_type;

		LMAT_ENSURE_INLINE
		map_multicol_reader(const Fun& fun, U, const Rd1& rd1, const Rd2& rd2)
		:m_fun(fun), m_rd1(rd1), m_rd2(rd2) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_fun, U(), m_rd1.col(j), m_rd2.col(j));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
		Rd2 m_rd2;
	};


	template<typename Fun, typename U, typename Rd1, typename Rd2, typename Rd3>
	class map_multicol_reader<Fun, U, Rd1, Rd2, Rd3> : public multicol_accessor_base
	{
		typedef typename Rd1::col_accessor_type col1_t;
		typedef typename Rd2::col_accessor_type col2_t;
		typedef typename Rd3::col_accessor_type col3_t;

	public:
		typedef map_vec_reader<Fun, U, col1_t, col2_t, col3_t> col_accessor_type;

		LMAT_ENSURE_INLINE
		map_multicol_reader(const Fun& fun, U, const Rd1& rd1, const Rd2& rd2, const Rd3& rd3)
		:m_fun(fun), m_rd1(rd1), m_rd2(rd2), m_rd3(rd3) { }

		LMAT_ENSURE_INLINE
		col_accessor_type col(index_t j) const
		{
			return col_accessor_type(m_fun, U(), m_rd1.col(j), m_rd2.col(j), m_rd3.col(j));
		}

	private:
		Fun m_fun;
		Rd1 m_rd1;
		Rd2 m_rd2;
		Rd3 m_rd3;
	};



}

#endif /* MAP_ACCESSORS_H_ */
