/**
 * @file mat_cast.h
 *
 * @brief Matrix value casting
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAT_CAST_H_
#define LIGHTMAT_MAT_CAST_H_

#include <light_mat/mateval/map_expr.h>

#define LMAT_DEFINE_MAT_CAST_FUN(T, Fun) \
	template<typename S, class SMat> \
	LMAT_ENSURE_INLINE \
	inline map_expr<cast_<T>, SMat> \
	Fun(const IEWiseMatrix<SMat, S>& smat) { \
		return make_map_expr(cast_<T>(), smat); }

namespace lmat
{
	/********************************************
	 *
	 *  definitions of functors and maps
	 *
	 ********************************************/

	template<typename T>
	struct cast_ { };

	template<typename S, typename T>
	struct fun_traits<cast_<T>, S>
	{
		typedef T result_type;
	};

	template<typename S, typename T>
	struct cast_fun
	{
		typedef T value_type;

		LMAT_ENSURE_INLINE
		T operator() (const S& s) const
		{
			return static_cast<T>(s);
		}
	};

	template<typename S, typename T>
	struct fun_map<cast_<T>, S>
	{
		typedef cast_fun<S, T> type;
	};


	/********************************************
	 *
	 *  matrix functions
	 *
	 ********************************************/

	template<typename S, class SMat, typename T>
	LMAT_ENSURE_INLINE
	inline map_expr<cast_<T>, SMat>
	cast(const IEWiseMatrix<SMat, S>& smat, type_<T>)
	{
		return make_map_expr(cast_<T>(), smat);
	}


	LMAT_DEFINE_MAT_CAST_FUN(double, to_f64)
	LMAT_DEFINE_MAT_CAST_FUN(float,  to_f32)

	LMAT_DEFINE_MAT_CAST_FUN( int8_t, to_i8)
	LMAT_DEFINE_MAT_CAST_FUN(uint8_t, to_u8)
	LMAT_DEFINE_MAT_CAST_FUN( int16_t, to_i16)
	LMAT_DEFINE_MAT_CAST_FUN(uint16_t, to_u16)
	LMAT_DEFINE_MAT_CAST_FUN( int32_t, to_i32)
	LMAT_DEFINE_MAT_CAST_FUN(uint32_t, to_u32)
	LMAT_DEFINE_MAT_CAST_FUN( int64_t, to_i64)
	LMAT_DEFINE_MAT_CAST_FUN(uint64_t, to_u64)

	LMAT_DEFINE_MAT_CAST_FUN( bool, to_bool)
	LMAT_DEFINE_MAT_CAST_FUN( mask_t<double>, to_f64m )
	LMAT_DEFINE_MAT_CAST_FUN( mask_t<float>,  to_f32m )

}

#endif /* MAT_CAST_H_ */
