/*
 * @file matlab_port.h
 *
 * @brief
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MATLAB_PORT_H_
#define LIGHTMAT_MATLAB_PORT_H_

#include <light_mat/matlab/marray.h>
#include <light_mat/matrix/matrix_classes.h>
#include <light_mat/matexpr/mat_arith.h>
#include <light_mat/common/arg_check.h>

// useful macros

#ifdef LMAT_DISABLE_MX_CHECK
#define LMAT_MX( idx, name, MTag, T ) auto name = lmat::matlab::as_(prhs[idx], lmat::matlab::MTag<T>());
#else
#define LMAT_MX( idx, name, MTag, T ) \
	const_marray _##name##_mx_(prhs[idx]); \
	if (lmat::matlab::MTag<T>::require_scalar) { \
		if (!_##name##_mx_.is_scalar()) { throw lmat::invalid_argument(#name " must be a scalar."); } \
	} \
	else { \
		if (!_##name##_mx_.is_matrix()) { throw lmat::invalid_argument(#name " must be a matrix."); } \
	} \
	if (!_##name##_mx_.is_of_type<T>()) { throw lmat::invalid_argument(#name " has invalid value type."); } \
	if (_##name##_mx_.is_sparse()) { throw lmat::invalid_argument(#name " can not be a sparse-matrix."); } \
	auto name = lmat::matlab::as_(prhs[idx], lmat::matlab::MTag<T>());
#endif

#define LMAT_MX_OUT( idx, name, MObj, MTag, T ) \
	plhs[idx] = MObj.get(); \
	auto name = lmat::matlab::as_(plhs[idx], lmat::matlab::MTag<T>());

#define LMAT_MX_NARGINCHK(L, H) \
		if (nrhs < L) throw lmat::invalid_argument("Not enough inputs."); \
		if (nrhs > H) throw lmat::invalid_argument("Too many inputs.");

#define LMAT_MX_NARGOUTCHK(L, H) \
		if (nlhs < L) throw lmat::invalid_argument("Not enough outputs."); \
		if (nlhs > H) throw lmat::invalid_argument("Too many outputs.");

void _lmat_mex_main_entry(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
template<typename T> void _lmat_mex_generic_main_entry(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

namespace lmat { namespace matlab {

	/********************************************
	 *
	 *  view mapping functions
	 *
	 ********************************************/

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_matrix<T> view2d(const_marray m)
	{
		return cref_matrix<T>(m.data<T>(), m.nrows(), m.ncolumns());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_matrix<T> view2d(marray m)
	{
		return ref_matrix<T>(m.data<T>(), m.nrows(), m.ncolumns());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_col<T> view_as_col(const_marray m)
	{
		return cref_col<T>(m.data<T>(), m.nelems());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_col<T> view_as_col(marray m)
	{
		return ref_col<T>(m.data<T>(), m.nelems());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_row<T> view_as_row(const_marray m)
	{
		return cref_row<T>(m.data<T>(), m.nelems());
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_row<T> view_as_row(marray m)
	{
		return ref_row<T>(m.data<T>(), m.nelems());
	}


	template<class Mat, typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const IMatrixXpr<Mat, T>& in)
	{
		marray r = marray::numeric_matrix<T>(in.nrows(), in.ncolumns());

		typedef ref_matrix<T, meta::nrows<Mat>::value, meta::ncols<Mat>::value> view_t;
		view_t view(r.data<T>(), in.nrows(), in.ncolumns());
		view = in;

		return r;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const T *src, index_t m, index_t n)
	{
		marray r = marray::numeric_matrix<T>(m, n);
		copy_vec(m * n, src, r.data<T>());
		return r;
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline marray to_marray(const double *src, index_t m, index_t n)
	{
		marray r = marray::double_matrix(m, n);
		copy_vec(m * n, src, r.data<T>());
		return r;
	}


	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	static typename std::enable_if<std::is_arithmetic<T>::value,
	marray>::type
	marray_like(const IMatrixXpr<Mat, T>& src)
	{
		return marray::numeric_matrix<T>(src.nrows(), src.ncolumns());
	}

	template<typename T, class Mat>
	LMAT_ENSURE_INLINE
	static marray marray_like(const IMatrixXpr<Mat, bool>& src)
	{
		return marray::logical_matrix(src.nrows(), src.ncolumns());
	}


	/********************************************
	 *
	 *  convenient tag based dispatching
	 *
	 ********************************************/

	template<typename T> struct mat_ { static const bool require_scalar = false; };
	template<typename T> struct col_ { static const bool require_scalar = false; };
	template<typename T> struct row_ { static const bool require_scalar = false; };
	template<typename T> struct sca_ { static const bool require_scalar = true; };

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_matrix<T> as_(const mxArray* mx, mat_<T>)
	{
		return view2d<T>(const_marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_matrix<T> as_(mxArray* mx, mat_<T>)
	{
		return view2d<T>(marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_col<T> as_(const mxArray* mx, col_<T>)
	{
		return view_as_col<T>(const_marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_col<T> as_(mxArray* mx, col_<T>)
	{
		return view_as_col<T>(marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline cref_row<T> as_(const mxArray* mx, row_<T>)
	{
		return view_as_row<T>(const_marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline ref_row<T> as_(mxArray* mx, row_<T>)
	{
		return view_as_row<T>(marray(mx));
	}

	template<typename T>
	LMAT_ENSURE_INLINE
	inline T as_(const mxArray* mx, sca_<T>)
	{
		return const_marray(mx).scalar<T>();
	}


	/********************************************
	 *
	 *  classes for dispatch
	 *
	 ********************************************/

	template<typename... AT>
	struct mex_main_dispatch;

#define _LMAT_MEX_DISPATCH_ENTRY(I) \
		case type_to_mclassid<AT##I>::value: \
			_lmat_mex_generic_main_entry<AT##I>(nlhs, plhs, nrhs, prhs); \
			return true

#define _LMAT_DEFINE_MEX_MAIN_DISPATCH(N) \
		template<LMAT_REPEAT_ARGS_P##N(typename AT)> \
		struct mex_main_dispatch<LMAT_REPEAT_ARGS_P##N(AT)> { \
			bool run(int idx, int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { \
				switch (mxGetClassID(prhs[idx])) { \
					LMAT_REPEAT_STATEMENTS_F##N(_LMAT_MEX_DISPATCH_ENTRY); \
					default: return false; \
				} \
			} };


	_LMAT_DEFINE_MEX_MAIN_DISPATCH(1)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(2)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(3)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(4)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(5)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(6)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(7)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(8)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(9)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(10)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(11)
	_LMAT_DEFINE_MEX_MAIN_DISPATCH(12)

} }

/********************************************
 *
 *  Macros for dispatching
 *
 ********************************************/

#define _LMAT_MEX_CATCH_PART(mname) \
		catch (lmat::invalid_argument& exc) { \
			mexErrMsgIdAndTxt( #mname ":invalidarg", exc.what() ); } \
		catch (lmat::invalid_operation& exc) { \
			mexErrMsgIdAndTxt( #mname ":invalidop", exc.what() ); } \
		catch (lmat::out_of_range& exc) { \
			mexErrMsgIdAndTxt( #mname ":outrange", exc.what() ); } \
		catch (std::exception& exc) { \
			mexErrMsgIdAndTxt( #mname ":stdexception", exc.what() ); }

#define LMAT_SIMPLE_MEX(mname) \
	void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { \
		try { \
			_lmat_mex_main_entry(nlhs, plhs, nrhs, prhs); \
		} \
		_LMAT_MEX_CATCH_PART(mname) } \
	void _lmat_mex_main_entry(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

#define LMAT_GENERIC_MEX(mname, idx, ...) \
	void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) { \
		try { \
			if (nrhs <= idx) throw lmat::invalid_argument("Not enough inputs"); \
			lmat::matlab::mex_main_dispatch<__VA_ARGS__> main_dispatcher; \
			if (!main_dispatcher.run(idx, nlhs, plhs, nrhs, prhs)) { \
				throw lmat::invalid_argument("Unsupported master value type"); \
			} \
		} \
		_LMAT_MEX_CATCH_PART(mname) } \
	template<typename T> \
	void _lmat_mex_generic_main_entry(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

#define LMAT_FP_MEX(mname, idx) LMAT_GENERIC_MEX(mname, idx, double, float)

#define LMAT_INT_MEX(mname, idx) \
	LMAT_GENERIC_MEX(mname, idx, int32_t, uint32_t, int64_t, uint64_t, int8_t, uint8_t, int16_t, uint16_t)

#define LMAT_NUM_MEX(mname, idx) \
	LMAT_GENERIC_MEX(mname, idx, double, float, int32_t, uint32_t, int64_t, uint64_t, int8_t, uint8_t, int16_t, uint16_t)

#define LMAT_FP_AND_BOOL_MEX(mname, idx) LMAT_GENERIC_MEX(mname, idx, double, float, bool)

#define LMAT_INT_AND_BOOL_MEX(mname, idx) \
	LMAT_GENERIC_MEX(mname, idx, int32_t, uint32_t, int64_t, uint64_t, int8_t, uint8_t, int16_t, uint16_t, bool)

#define LMAT_NUM_AND_BOOL_MEX(mname, idx) \
	LMAT_GENERIC_MEX(mname, idx, double, float, int32_t, uint32_t, int64_t, uint64_t, int8_t, uint8_t, int16_t, uint16_t, bool)

#endif /* MATLAB_PORT_H_ */
