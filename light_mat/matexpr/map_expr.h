/**
 * @file map_expr.h
 *
 * Map expression classes
 * 
 * @author Dahua Lin 
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MAP_EXPR_H_
#define LIGHTMAT_MAP_EXPR_H_

#include <light_mat/matexpr/matexpr_fwd.h>
#include <light_mat/math/fun_tags.h>
#include <light_mat/mateval/ewise_eval.h>

#include "internal/map_expr_internal.h"

namespace lmat
{

	/********************************************
	 *
	 *  matrix traits
	 *
	 ********************************************/

	template<typename Tag, typename Arg1>
	struct matrix_traits<map_expr<Tag, Arg1> >
	{
		typedef typename internal::mapexpr_helper_map<Arg1>::type helper;

		static const int num_dimensions = 2;
		static const int ct_num_rows = helper::ct_rows;
		static const int ct_num_cols = helper::ct_cols;
		static const bool is_readonly = true;

		typedef typename helper::vtype1 arg1_value_type;

		typedef typename fun_traits<Tag, arg1_value_type>::result_type value_type;
		typedef typename helper::domain domain;
	};

	template<typename Tag, typename Arg1, typename Arg2>
	struct matrix_traits<map_expr<Tag, Arg1, Arg2> >
	{
		typedef typename internal::mapexpr_helper_map<Arg1, Arg2>::type helper;

		static const int num_dimensions = 2;
		static const int ct_num_rows = helper::ct_rows;
		static const int ct_num_cols = helper::ct_cols;
		static const bool is_readonly = true;

		typedef typename helper::vtype1 arg1_value_type;
		typedef typename helper::vtype2 arg2_value_type;

		typedef typename fun_traits<Tag,
				arg1_value_type, arg2_value_type>::result_type value_type;
		typedef typename helper::domain domain;
	};

	template<typename Tag, typename Arg1, typename Arg2, typename Arg3>
	struct matrix_traits<map_expr<Tag, Arg1, Arg2, Arg3> >
	{
		typedef typename internal::mapexpr_helper_map<Arg1, Arg2, Arg3>::type helper;

		static const int num_dimensions = 2;
		static const int ct_num_rows = helper::ct_rows;
		static const int ct_num_cols = helper::ct_cols;
		static const bool is_readonly = true;

		typedef typename helper::vtype1 arg1_value_type;
		typedef typename helper::vtype2 arg2_value_type;
		typedef typename helper::vtype3 arg3_value_type;

		typedef typename fun_traits<Tag,
				arg1_value_type, arg2_value_type, arg3_value_type>::result_type value_type;
		typedef typename helper::domain domain;
	};



	/********************************************
	 *
	 *  expression classes
	 *
	 ********************************************/

	template<typename Tag, typename Arg1>
	class map_expr<Tag, Arg1>
	: public IMatrixXpr<map_expr<Tag, Arg1>,
	  typename matrix_traits<map_expr<Tag, Arg1> >::value_type>
	{
		typedef typename internal::mapexpr_helper_map<Arg1>::type helper;

	public:
		typedef Tag tag_type;
		typedef Arg1 arg1_type;
		typedef typename helper::value_type value_type;
		typedef typename helper::shape_type shape_type;

		LMAT_ENSURE_INLINE
		map_expr(const Tag& tag, const Arg1& a1)
		: m_tag(tag), m_shape(helper::get_shape(a1))
		, m_arg1(a1) { }

		LMAT_ENSURE_INLINE
		const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE
		const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		index_t shape() const
		{
			return m_shape;
		}

	private:
		Tag m_tag;
		shape_type m_shape;
		const Arg1& m_arg1;
	};


	template<typename Tag, typename Arg1, typename Arg2>
	class map_expr<Tag, Arg1, Arg2>
	: public IMatrixXpr<map_expr<Tag, Arg1, Arg2>,
	  typename matrix_traits<map_expr<Tag, Arg1, Arg2> >::value_type>
	{
		typedef typename internal::mapexpr_helper_map<Arg1, Arg2>::type helper;

	public:
		typedef Tag tag_type;
		typedef Arg1 arg1_type;
		typedef Arg2 arg2_type;
		typedef typename helper::value_type value_type;
		typedef typename helper::shape_type shape_type;

		LMAT_ENSURE_INLINE
		map_expr(const Tag& tag, const Arg1& a1, const Arg2& a2)
		: m_tag(tag), m_shape(helper::get_shape(a1, a2))
		, m_arg1(a1), m_arg2(a2) { }

		LMAT_ENSURE_INLINE
		const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE
		const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE
		const Arg2& arg2() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		index_t shape() const
		{
			return m_shape;
		}

	private:
		Tag m_tag;
		shape_type m_shape;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};


	template<typename Tag, typename Arg1, typename Arg2, typename Arg3>
	class map_expr<Tag, Arg1, Arg2, Arg3>
	: public IMatrixXpr<map_expr<Tag, Arg1, Arg2, Arg3>,
	  typename matrix_traits<map_expr<Tag, Arg1, Arg2, Arg3> >::value_type>
	{
		typedef typename internal::mapexpr_helper_map<Arg1, Arg2, Arg3>::type helper;

	public:
		typedef Tag tag_type;
		typedef Arg1 arg1_type;
		typedef Arg2 arg2_type;
		typedef Arg3 arg3_type;
		typedef typename helper::value_type value_type;
		typedef typename helper::shape_type shape_type;

		LMAT_ENSURE_INLINE
		map_expr(const Tag& tag, const Arg1& a1, const Arg2& a2, const Arg3& a3)
		: m_tag(tag), m_shape(helper::get_shape(a1, a2, a3))
		, m_arg1(a1), m_arg2(a2), m_arg3(a3) { }

		LMAT_ENSURE_INLINE
		const Tag& tag() const
		{
			return m_tag;
		}

		LMAT_ENSURE_INLINE
		const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE
		const Arg2& arg2() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE
		const Arg3& arg3() const
		{
			return m_arg3;
		}

		LMAT_ENSURE_INLINE
		index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t ncolumns() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE
		index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE
		index_t shape() const
		{
			return m_shape;
		}

	private:
		Tag m_tag;
		shape_type m_shape;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
		const Arg3& m_arg3;
	};


	/********************************************
	 *
	 *  readers
	 *
	 ********************************************/


	/********************************************
	 *
	 *  evaluation
	 *
	 ********************************************/

	template<typename Tag, typename... Args, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const map_expr<Tag, Args...>& expr,
			IRegularMatrix<DMat, typename matrix_traits<map_expr<Tag, Args...> >::value_type>& dst)
	{

	}

}

#endif 
