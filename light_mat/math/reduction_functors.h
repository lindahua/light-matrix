/**
 * @file reduction_functors.h
 *
 * Basic reduction functors
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_REDUCTION_FUNCTORS_H_
#define LIGHTMAT_REDUCTION_FUNCTORS_H_

#include <light_mat/common/arg_check.h>
#include <light_mat/math/functor_base.h>
#include <light_mat/math/math_base.h>


/************************************************
 *
 *  Macros to declare reduction tags
 *
 ************************************************/

// reduction on generic types

#define LMAT_DECLARE_GENERIC_UNARY_REDUCTION( Name, InterT, ResultT ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 1> { static const bool value = false; }; \
	template<typename T> \
	struct reduction_result<Name##_t, T> { \
		typedef InterT intermediate_type; \
		typedef ResultT type; };

#define LMAT_DECLARE_GENERIC_BINARY_REDUCTION( Name, InterT, ResultT ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 2> { static const bool value = false; }; \
	template<typename T> \
	struct reduction_result<Name##_t, T, T> { \
		typedef InterT intermediate_type; \
		typedef ResultT type; };

#define LMAT_DECLARE_GENERIC_SIMPLE_UNARY_REDUCTION( Name ) \
	LMAT_DECLARE_GENERIC_UNARY_REDUCTION( Name, T, T )

#define LMAT_DECLARE_GENERIC_SIMPLE_BINARY_REDUCTION( Name ) \
	LMAT_DECLARE_GENERIC_BINARY_REDUCTION( Name, T, T )


#define LMAT_DECLARE_REDUCTION_WITHPOST( Name ) \
	template<> struct reduction_with_post<Name##_t> { static const bool value = true; };


// reduction on real numbers

#define LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( Name ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 1> { static const bool value = false; }; \
	template<> struct reduction_result<Name##_t, float> { \
		typedef float intermediate_type; \
		typedef float type; }; \
	template<> struct reduction_result<Name##_t, double> { \
		typedef double intermediate_type; \
		typedef double type; };

#define LMAT_DECLARE_REAL_SIMPLE_BINARY_REDUCTION( Name ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 2> { static const bool value = false; }; \
	template<> struct reduction_result<Name##_t, float, float> { \
		typedef float intermediate_type; \
		typedef float type; }; \
	template<> struct reduction_result<Name##_t, double, double> { \
		typedef double intermediate_type; \
		typedef double type; };

#define LMAT_DECLARE_REAL_MEDIATED_UNARY_REDUCTION( Name, InterTempl ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 1> { static const bool value = false; }; \
	template<> struct reduction_result<Name##_t, float> { \
		typedef InterTempl<float> intermediate_type; \
		typedef float type; }; \
	template<> struct reduction_result<Name##_t, double> { \
		typedef InterTempl<float> intermediate_type; \
		typedef double type; };

#define LMAT_DECLARE_REAL_MEDIATED_BINARY_REDUCTION( Name, InterTempl ) \
	struct Name##_t { }; \
	template<> struct is_reduction_tag<Name##_t, 2> { static const bool value = false; }; \
	template<> struct reduction_result<Name##_t, float, float> { \
		typedef InterTempl<float> intermediate_type; \
		typedef float type; }; \
	template<> struct reduction_result<Name##_t, double, double> { \
		typedef InterTempl<double> intermediate_type; \
		typedef double type; };


/************************************************
 *
 *  Macros to define reduction functors
 *
 ************************************************/

#define LMAT_DEFINE_UNARY_REDUCTOR(Name, EmptyVal, TransformExpr, CombExpr, GetExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef typename reduction_result<Name##_t, T>::intermediate_type intermediate_type; \
		typedef typename reduction_result<Name##_t, T>::type result_type; \
		LMAT_ENSURE_INLINE \
		Name##_fun() { } \
		LMAT_ENSURE_INLINE \
		Name##_fun( Name##_t ) { } \
		LMAT_ENSURE_INLINE \
		T empty_value() const { return EmptyVal; } \
		LMAT_ENSURE_INLINE \
		T transform(const T& x) const { return TransformExpr; } \
		LMAT_ENSURE_INLINE \
		T combine(const intermediate_type& a, const intermediate_type& b) const { return CombExpr; } \
		LMAT_ENSURE_INLINE \
		result_type get(const intermediate_type& a, const index_t& n) const { return GetExpr; } \
	}; \
	template<typename T> \
	struct reduction_fun<Name##_t, scalar_ker, T> { \
		typedef Name##_fun<T> type; \
	};

#define LMAT_DEFINE_BINARY_REDUCTOR(Name, EmptyVal, TransformExpr, CombExpr, GetExpr ) \
	template<typename T> \
	struct Name##_fun { \
		typedef typename reduction_result<Name##_t, T, T>::intermediate_type intermediate_type; \
		typedef typename reduction_result<Name##_t, T, T>::type result_type; \
		LMAT_ENSURE_INLINE \
		Name##_fun() { } \
		LMAT_ENSURE_INLINE \
		Name##_fun( Name##_t ) { } \
		LMAT_ENSURE_INLINE \
		result_type empty_value() const { return EmptyVal; } \
		LMAT_ENSURE_INLINE \
		intermediate_type transform(const T& x, const T& y) const { return TransformExpr; } \
		LMAT_ENSURE_INLINE \
		intermediate_type combine(const intermediate_type& a, const intermediate_type& b) const { return CombExpr; } \
		LMAT_ENSURE_INLINE \
		result_type get(const intermediate_type& a, const index_t& n) const { return GetExpr; } \
	}; \
	template<typename T> \
	struct reduction_fun<Name##_t, scalar_ker, T, T> { \
		typedef Name##_fun<T> type; \
	};


namespace lmat
{

	/********************************************
	 *
	 *  Basic devices
	 *
	 ********************************************/

	template<typename Tag, int NArgs>
	struct is_reduction_tag
	{
		static const bool value = false;
	};

	template<typename Tag>
	struct reduction_with_post
	{
		static const bool value = false;
	};

	template<typename Op, typename T1, typename T2=nil_t, typename T3=nil_t>
	struct reduction_result;

	template<typename Op, typename Ker, typename T1, typename T2=nil_t, typename T3=nil_t>
	struct reduction_fun;


	/********************************************
	 *
	 *  reduction tags
	 *
	 ********************************************/

	namespace internal
	{
		template<typename T> struct nrmdot_media;
	}


	LMAT_DECLARE_GENERIC_SIMPLE_UNARY_REDUCTION( sum )
	LMAT_DECLARE_GENERIC_SIMPLE_UNARY_REDUCTION( maximum )
	LMAT_DECLARE_GENERIC_SIMPLE_UNARY_REDUCTION( minimum )

	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( mean )
	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( L1norm )
	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( sqL2norm )
	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( L2norm )
	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( Linfnorm )

	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( logsum )
	LMAT_DECLARE_REAL_SIMPLE_UNARY_REDUCTION( entropy )

	LMAT_DECLARE_REAL_SIMPLE_BINARY_REDUCTION( dot )
	LMAT_DECLARE_REAL_MEDIATED_BINARY_REDUCTION( nrmdot, internal::nrmdot_media )

	// specify the ones that require post processing

	LMAT_DECLARE_REDUCTION_WITHPOST( mean )
	LMAT_DECLARE_REDUCTION_WITHPOST( L2norm )
	LMAT_DECLARE_REDUCTION_WITHPOST( entropy )
	LMAT_DECLARE_REDUCTION_WITHPOST( nrmdot )


	/********************************************
	 *
	 *  reduction functors
	 *
	 ********************************************/

	template<typename T>
	inline T no_empty_value(const char *msg)
	{
		throw invalid_operation(msg);
	}

	// sum, maximum, minimum, and mean

	LMAT_DEFINE_UNARY_REDUCTOR( sum,
			T(0),
			x,
			a + b,
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( maximum,
			no_empty_value<T>("maximum can not be applied to empty arrays."),
			x,
			(math::max)(a, b),
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( minimum,
			no_empty_value<T>("minimum can not be applied to empty arrays."),
			x,
			(math::min)(a, b),
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( mean,
			no_empty_value<T>("mean can not be applied to empty arrays."),
			x,
			a + b,
			a / static_cast<T>(n))

	// norms

	LMAT_DEFINE_UNARY_REDUCTOR( L1norm,
			T(0),
			math::abs(x),
			a + b,
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( sqL2norm,
			T(0),
			math::sqr(x),
			a + b,
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( L2norm,
			T(0),
			math::sqr(x),
			a + b,
			math::sqrt(a))

	LMAT_DEFINE_UNARY_REDUCTOR( Linfnorm,
			T(0),
			math::abs(x),
			(math::max)(a, b),
			a)


	// logsum & entropy

	namespace math
	{
		LMAT_ENSURE_INLINE inline float xlogx(float x) { return x > 0 ? x * math::log(x) : 0.f; }
		LMAT_ENSURE_INLINE inline double xlogx(double x) { return x > 0 ? x * math::log(x) : 0.0; }
	}

	LMAT_DEFINE_UNARY_REDUCTOR( logsum,
			T(0),
			math::log(x),
			a + b,
			a)

	LMAT_DEFINE_UNARY_REDUCTOR( entropy,
			T(0),
			math::xlogx(x),
			a + b,
			- a)

	// dot & nrmdot

	LMAT_DEFINE_BINARY_REDUCTOR( dot,
			T(0),
			x * y,
			a + b,
			a)


	namespace internal
	{
		template<typename T>
		struct nrmdot_media
		{
			T sum_xx;
			T sum_yy;
			T sum_xy;

			LMAT_ENSURE_INLINE
			T get() const
			{
				return sum_xy / (math::sqrt(sum_xx) * math::sqrt(sum_yy));
			}

			LMAT_ENSURE_INLINE
			nrmdot_media operator + (const nrmdot_media& b) const
			{
				nrmdot_media r;
				r.sum_xx = sum_xx + b.sum_xx;
				r.sum_yy = sum_yy + b.sum_yy;
				r.sum_xy = sum_xy + b.sum_xy;
				return r;
			}

			LMAT_ENSURE_INLINE
			static nrmdot_media from(const T& x, const T& y)
			{
				nrmdot_media r;
				r.sum_xx = x * x;
				r.sum_yy = y * y;
				r.sum_xy = x * y;
				return r;
			}
		};
	}

	LMAT_DEFINE_BINARY_REDUCTOR( nrmdot,
			T(0),
			internal::nrmdot_media<T>::from(x, y),
			a + b,
			a.get())

}

#endif





