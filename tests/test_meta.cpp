/**
 * @file test_meta.cpp
 *
 * @brief Test meta programming constructs
 *
 * @author Dahua Lin
 */


#include "test_base.h"

#include <light_mat/common/meta_base.h>

using namespace lmat;

SIMPLE_CASE( meta_calc, logical )
{
	ASSERT_CT_VALUE( true_type, true );
	ASSERT_CT_VALUE( false_type, false );

	typedef meta::and_<true_type,  true_type>  and_tt;
	typedef meta::and_<true_type,  false_type> and_tf;
	typedef meta::and_<false_type, true_type>  and_ft;
	typedef meta::and_<false_type, false_type> and_ff;

	ASSERT_CT_VALUE( and_tt, true );
	ASSERT_CT_VALUE( and_tf, false );
	ASSERT_CT_VALUE( and_ft, false );
	ASSERT_CT_VALUE( and_ff, false );

	typedef meta::or_<true_type,  true_type>  or_tt;
	typedef meta::or_<true_type,  false_type> or_tf;
	typedef meta::or_<false_type, true_type>  or_ft;
	typedef meta::or_<false_type, false_type> or_ff;

	ASSERT_CT_VALUE( or_tt, true );
	ASSERT_CT_VALUE( or_tf, true );
	ASSERT_CT_VALUE( or_ft, true );
	ASSERT_CT_VALUE( or_ff, false );
}

SIMPLE_CASE( meta_calc, arith )
{
	typedef fix_int<2> i2;
	typedef fix_int<3> i3;

	typedef meta::add_<i2, i3> add_r;
	typedef meta::sub_<i2, i3> sub_r;
	typedef meta::mul_<i2, i3> mul_r;
	typedef meta::max_<i2, i3> max_r;
	typedef meta::min_<i2, i3> min_r;

	ASSERT_CT_VALUE( add_r, 5 );
	ASSERT_CT_VALUE( sub_r, -1 );
	ASSERT_CT_VALUE( mul_r, 6 );
	ASSERT_CT_VALUE( max_r, 3 );
	ASSERT_CT_VALUE( min_r, 2 );
}


SIMPLE_CASE( typelist, N1 )
{
	typedef meta::type_list<int> list_t;

	ASSERT_CT_VALUE ( meta::len_<list_t>, 1 );

	typedef typename meta::get_<list_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, int );
}


SIMPLE_CASE( typelist, N2 )
{
	typedef meta::type_list<int, char> list_t;

	ASSERT_CT_VALUE ( meta::len_<list_t>, 2 );

	typedef typename meta::get_<list_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, int );

	typedef typename meta::get_<list_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, char );
}


SIMPLE_CASE( typelist, N3 )
{
	typedef meta::type_list<int, char, float> list_t;

	ASSERT_CT_VALUE ( meta::len_<list_t>, 3 );

	typedef typename meta::get_<list_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, int );

	typedef typename meta::get_<list_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, char );

	typedef typename meta::get_<list_t, 2>::type T2;
	ASSERT_SAME_TYPE ( T2, float );
}


SIMPLE_CASE( typelist, N4 )
{
	typedef meta::type_list<int, char, float, double> list_t;

	ASSERT_CT_VALUE ( meta::len_<list_t>, 4 );

	typedef typename meta::get_<list_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, int );

	typedef typename meta::get_<list_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, char );

	typedef typename meta::get_<list_t, 2>::type T2;
	ASSERT_SAME_TYPE ( T2, float );

	typedef typename meta::get_<list_t, 3>::type T3;
	ASSERT_SAME_TYPE ( T3, double );
}


template<typename T>
struct my_map
{
	typedef type_<T> type;
};


SIMPLE_CASE( typemap, N1 )
{
	typedef meta::type_list<int> list_t;
	typedef meta::map_<my_map, list_t>::type result_t;

	typedef typename meta::get_<result_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, type_<int> );
}

SIMPLE_CASE( typemap, N2 )
{
	typedef meta::type_list<int, char> list_t;
	typedef meta::map_<my_map, list_t>::type result_t;

	typedef typename meta::get_<result_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, type_<int> );

	typedef typename meta::get_<result_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, type_<char> );
}

SIMPLE_CASE( typemap, N3 )
{
	typedef meta::type_list<int, char, float> list_t;
	typedef meta::map_<my_map, list_t>::type result_t;

	typedef typename meta::get_<result_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, type_<int> );

	typedef typename meta::get_<result_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, type_<char> );

	typedef typename meta::get_<result_t, 2>::type T2;
	ASSERT_SAME_TYPE ( T2, type_<float> );
}


SIMPLE_CASE( typemap, N4 )
{
	typedef meta::type_list<int, char, float, double> list_t;
	typedef meta::map_<my_map, list_t>::type result_t;

	typedef typename meta::get_<result_t, 0>::type T0;
	ASSERT_SAME_TYPE ( T0, type_<int> );

	typedef typename meta::get_<result_t, 1>::type T1;
	ASSERT_SAME_TYPE ( T1, type_<char> );

	typedef typename meta::get_<result_t, 2>::type T2;
	ASSERT_SAME_TYPE ( T2, type_<float> );

	typedef typename meta::get_<result_t, 3>::type T3;
	ASSERT_SAME_TYPE ( T3, type_<double> );
}



BEGIN_TPACK( meta_calc )
	ADD_SIMPLE_CASE( meta_calc, logical )
	ADD_SIMPLE_CASE( meta_calc, arith )
END_TPACK

BEGIN_TPACK( typelists )
	ADD_SIMPLE_CASE( typelist, N1 )
	ADD_SIMPLE_CASE( typelist, N2 )
	ADD_SIMPLE_CASE( typelist, N3 )
	ADD_SIMPLE_CASE( typelist, N4 )
END_TPACK

BEGIN_TPACK( typemaps )
	ADD_SIMPLE_CASE( typemap, N1 )
	ADD_SIMPLE_CASE( typemap, N2 )
	ADD_SIMPLE_CASE( typemap, N3 )
	ADD_SIMPLE_CASE( typemap, N4 )
END_TPACK



BEGIN_MAIN_SUITE
	ADD_TPACK( meta_calc )
	ADD_TPACK( typelists )
	ADD_TPACK( typemaps )
END_MAIN_SUITE



