/**
 * @file common_shape_infer.h
 *
 * Inference of common shape for map expression
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAP_EXPR_INTERNAL_H_
#define LIGHTMAT_MAP_EXPR_INTERNAL_H_

#include <light_mat/matexpr/matexpr_fwd.h>

#include <light_mat/matrix/matrix_properties.h>
#include <light_mat/mateval/multicol_accessors.h>
#include <light_mat/mateval/fun_maps.h>

namespace lmat { namespace internal {

	template<class Mat> struct _mat_arg { };
	template<typename T> struct _sca_arg { };

	template<typename A>
	struct _to_qarg
	{
		typedef typename meta::if_<meta::is_mat_xpr<A>,
				_mat_arg<A>,
				_sca_arg<A> >::type type;
	};

	template<typename... QArg> struct mapexpr_helper;

	/********************************************
	 *
	 *  mapexpr helper map
	 *
	 ********************************************/

	template<typename... Arg>
	struct mapexpr_helper_map;

	template<typename A1>
	struct mapexpr_helper_map<A1>
	{
		typedef mapexpr_helper<
					typename _to_qarg<A1>::value_type
				> type;
	};

	template<typename A1, typename A2>
	struct mapexpr_helper_map<A1, A2>
	{
		typedef mapexpr_helper<
					typename _to_qarg<A1>::value_type,
					typename _to_qarg<A2>::value_type
				> type;
	};

	template<typename A1, typename A2, typename A3>
	struct mapexpr_helper_map<A1, A2, A3>
	{
		typedef mapexpr_helper<
					typename _to_qarg<A1>::value_type,
					typename _to_qarg<A2>::value_type,
					typename _to_qarg<A3>::value_type
				> type;
	};


	/********************************************
	 *
	 *  mapexpr helper classes
	 *
	 ********************************************/

	// unary

	template<typename A1>
	struct mapexpr_helper<_mat_arg<A1> >
	{
		typedef typename meta::domain_of<A1>::type domain;

		typedef typename matrix_traits<A1>::value_type vtype1;
		static const int ct_nrows = meta::nrows<A1>::value;
		static const int ct_ncols = meta::ncols<A1>::value;

		typedef matrix_shape<ct_nrows, ct_ncols> shape_type;

		LMAT_ENSURE_INLINE
		static shape_type get_shape(const A1& a1)
		{
			return a1.shape();
		}
	};


	// binary

	template<typename A1, typename A2>
	struct mapexpr_helper<_mat_arg<A1>, _mat_arg<A2> >
	{
		typedef typename meta::common_domain<A1, A2>::type domain;

		typedef typename matrix_traits<A1>::value_type vtype1;
		typedef typename matrix_traits<A2>::value_type vtype2;

		static const int ct_nrows = meta::common_nrows<A1, A2>::value;
		static const int ct_ncols = meta::common_ncols<A1, A2>::value;

		typedef matrix_shape<ct_nrows, ct_ncols> shape_type;

		LMAT_ENSURE_INLINE
		static shape_type get_shape(const A1& a1, const A2& a2)
		{
			return common_shape(a1, a2);
		}
	};

	template<typename A1, typename A2>
	struct mapexpr_helper<_mat_arg<A1>, _sca_arg<A2> >
	{
		typedef typename meta::domain_of<A1>::type domain;

		typedef typename matrix_traits<A1>::value_type vtype1;
		typedef A2 vtype2;

		static const int ct_nrows = meta::nrows<A1>::value;
		static const int ct_ncols = meta::ncols<A1>::value;

		typedef matrix_shape<ct_nrows, ct_ncols> shape_type;

		LMAT_ENSURE_INLINE
		static shape_type get_shape(const A1& a1, const A2& a2)
		{
			return a1.shape();
		}
	};

	template<typename A1, typename A2>
	struct mapexpr_helper<_sca_arg<A1>, _mat_arg<A2> >
	{
		typedef typename meta::domain_of<A2>::type domain;

		typedef A1 vtype1;
		typedef typename matrix_traits<A2>::value_type vtype2;

		static const int ct_nrows = meta::nrows<A2>::value;
		static const int ct_ncols = meta::ncols<A2>::value;

		typedef matrix_shape<ct_nrows, ct_ncols> shape_type;

		LMAT_ENSURE_INLINE
		static shape_type get_shape(const A1& a1, const A2& a2)
		{
			return a2.shape();
		}
	};


	// ternary

	template<typename A1, typename A2, typename A3>
	struct mapexpr_helper<_mat_arg<A1>, _mat_arg<A2>, _mat_arg<A3> >
	{
		typedef typename meta::common_domain<A1, A2, A3>::type domain;

		typedef typename matrix_traits<A1>::value_type vtype1;
		typedef typename matrix_traits<A2>::value_type vtype2;
		typedef typename matrix_traits<A3>::value_type vtype3;

		static const int ct_nrows = meta::common_nrows<A1, A2, A3>::value;
		static const int ct_ncols = meta::common_ncols<A1, A2, A3>::value;

		typedef matrix_shape<ct_nrows, ct_ncols> shape_type;

		LMAT_ENSURE_INLINE
		static shape_type get_shape(const A1& a1, const A2& a2, const A3& a3)
		{
			return common_shape(a1, a2, a3);
		}
	};


	/********************************************
	 *
	 *  vec reader classes
	 *
	 ********************************************/

	template<typename Tag, typename U, typename... Args> class map_vec_reader;
	template<typename Tag, typename U, typename... Args> class map_vec_reader;

	template<typename QArg, typename U> struct q_vec_reader_map;
	template<typename QArg, typename U> struct q_multicol_reader_map;


#define LMAT_MAP_READER_DEFS_S1( RdType ) \
		typedef typename _to_qarg<A1>::type qA1; \
		typedef mapexpr_helper<qA1> helper_t; \
		typedef q_##RdType##_reader_map<qA1, atag> map1_t; \
		typedef typename fun_traits<Tag, typename helper_t::vtype1>::result_type scalar_t; \
		typedef typename fun_map<Tag, typename helper_t::vtype1>::type fun_t;

#define LMAT_MAP_READER_DEFS_P1( RdType ) \
		typedef typename _to_qarg<A1>::type qA1; \
		typedef mapexpr_helper<qA1> helper_t; \
		typedef q_##RdType##_reader_map<qA1, atag> map1_t; \
		typedef typename fun_traits<Tag, typename helper_t::vtype1>::result_type scalar_t; \
		typedef math::simd_pack<scalar_t, SKind> pack_t; \
		typedef typename fun_map<Tag, typename helper_t::vtype1>::type fun_t;

#define LMAT_MAP_READER_DEFS_S2( RdType ) \
		typedef typename _to_qarg<A1>::type qA1; \
		typedef typename _to_qarg<A2>::type qA2; \
		typedef mapexpr_helper<qA1, qA2> helper_t; \
		typedef q_##RdType##_reader_map<qA1, atag> map1_t; \
		typedef q_##RdType##_reader_map<qA2, atag> map2_t; \
		typedef typename fun_traits<Tag, typename helper_t::vtype1, typename helper_t::vtype2>::result_type scalar_t; \
		typedef typename fun_map<Tag, typename helper_t::vtype1, typename helper_t::vtype2>::type fun_t;

#define LMAT_MAP_READER_DEFS_P2( RdType ) \
		typedef typename _to_qarg<A1>::type qA1; \
		typedef typename _to_qarg<A2>::type qA2; \
		typedef mapexpr_helper<qA1, qA2> helper_t; \
		typedef q_##RdType##_reader_map<qA1, atag> map1_t; \
		typedef q_##RdType##_reader_map<qA2, atag> map2_t; \
		typedef typename fun_traits<Tag, typename helper_t::vtype1, typename helper_t::vtype2>::result_type scalar_t; \
		typedef math::simd_pack<scalar_t, SKind> pack_t; \
		typedef typename fun_map<Tag, typename helper_t::vtype1, typename helper_t::vtype2>::type fun_t;




	template<typename Tag, typename A1>
	class map_vec_reader<Tag, atags::scalar, A1> : public scalar_vec_accessor_base
	{
		typedef atags::scalar atag;
		LMAT_MAP_READER_DEFS_S1(vec)

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(Tag, const A1& a1)
		: m_rd1(map1_t::get(a1)) { }

		LMAT_ENSURE_INLINE
		scalar_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i));
		}

	private:
		fun_t m_fun;
		typename map1_t::type m_rd1;
	};


	template<typename Tag, typename T, typename SKind, typename A1>
	class map_vec_reader<Tag, atags::simd<T, SKind>, A1> : public simd_vec_accessor_base
	{
		typedef atags::simd<T, SKind> atag;
		LMAT_MAP_READER_DEFS_P1(vec)

	public:
		LMAT_ENSURE_INLINE
		map_vec_reader(Tag, const A1& a1)
		: m_rd1(map1_t::get(a1)) { }

		LMAT_ENSURE_INLINE
		scalar_t scalar(index_t i) const
		{
			return m_fun(m_rd1.scalar(i));
		}

		LMAT_ENSURE_INLINE
		pack_t pack(index_t i) const
		{
			return m_fun(m_rd1.pack(i));
		}

	private:
		fun_t m_fun;
		typename map1_t::type m_rd1;
	};


	template<typename T, typename U>
	struct q_vec_reader_map<_sca_arg<T>, U>
	{
		typedef single_reader<T, U> type;

		LMAT_ENSURE_INLINE
		static type get(const T& a)
		{
			return single_reader<T, U>(a);
		}
	};

	template<typename A, typename U>
	struct q_vec_reader_map<_mat_arg<A>, U>
	{
		typedef typename internal::vec_reader_map<A, U>::type type;

		LMAT_ENSURE_INLINE
		static type get(const A& a)
		{
			return internal::vec_reader_map<A, U>::get(a);
		}
	};

	template<typename Tag, typename A1, typename U>
	struct q_vec_reader_map<_mat_arg<map_expr<Tag, A1> >, U>
	{
		typedef map_expr<Tag, A1> expr_t;
		typedef typename internal::vec_reader_map<expr_t, U>::type type;

		LMAT_ENSURE_INLINE
		static type get(const expr_t& a)
		{
			return map_vec_reader<Tag, U, A1>(a.tag(), a.arg1());
		}
	};

} }

#endif 
