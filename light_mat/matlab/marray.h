/*
 * @file marray.h
 *
 * @brief A light-weight wrapper of MATLAB array
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_MARRAY_H_
#define LIGHTMAT_MARRAY_H_

#include <mex.h>

#include <light_mat/common/basic_defs.h>

namespace lmat { namespace matlab {


	/********************************************
	 *
	 *  type maps
	 *
	 ********************************************/

	template<mxClassID C> struct mclassid_to_type;
	template<typename T> struct type_to_mclassid;

#define LMAT_DEFINE_MATLAB_TYPEMAP(C, T) \
	template<> struct mclassid_to_type<C> { typedef T type; }; \
	template<> struct type_to_mclassid<T> { static const mxClassID value = C; };

	LMAT_DEFINE_MATLAB_TYPEMAP(mxDOUBLE_CLASS, double)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxSINGLE_CLASS, float)

	LMAT_DEFINE_MATLAB_TYPEMAP(mxINT64_CLASS,  int64_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxUINT64_CLASS, uint64_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxINT32_CLASS,  int32_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxUINT32_CLASS, uint32_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxINT16_CLASS,  int16_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxUINT16_CLASS, uint16_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxINT8_CLASS,   int8_t)
	LMAT_DEFINE_MATLAB_TYPEMAP(mxUINT8_CLASS,  uint8_t)

	LMAT_DEFINE_MATLAB_TYPEMAP(mxLOGICAL_CLASS, bool)


	/********************************************
	 *
	 *  const_marray & marray classes
	 *
	 ********************************************/

	class const_marray
	{
	public:
		LMAT_ENSURE_INLINE
		const_marray(const mxArray *pm)
		: m_pm(const_cast<mxArray*>(pm))
		{ }

		LMAT_ENSURE_INLINE const mxArray *get() const
		{
			return m_pm;
		}

	public:
		// Type checking

		LMAT_ENSURE_INLINE mxClassID class_id() const
		{
			return mxGetClassID(m_pm);
		}

		LMAT_ENSURE_INLINE const char* class_name() const
		{
			return mxGetClassName(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_double() const
		{
			return mxIsDouble(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_single() const
		{
			return mxIsSingle(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_numeric() const
		{
			return mxIsNumeric(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_int64() const
		{
			return mxIsInt64(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_uint64() const
		{
			return mxIsUint64(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_int32() const
		{
			return mxIsInt32(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_uint32() const
		{
			return mxIsUint32(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_int16() const
		{
			return mxIsInt16(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_uint16() const
		{
			return mxIsUint16(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_int8() const
		{
			return mxIsInt8(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_uint8() const
		{
			return mxIsUint8(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_char() const
		{
			return mxIsChar(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_logical() const
		{
			return mxIsLogical(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_logical_scalar() const
		{
			return mxIsLogicalScalar(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_struct() const
		{
			return mxIsStruct(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_cell() const
		{
			return mxIsCell(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_empty() const
		{
			return mxIsEmpty(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_sparse() const
		{
			return mxIsSparse(m_pm);
		}

		LMAT_ENSURE_INLINE bool is_complex() const
		{
			return mxIsComplex(m_pm);
		}

	public:
		// Size

		LMAT_ENSURE_INLINE int ndims() const
		{
			return (int)mxGetNumberOfDimensions(m_pm);
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return (index_t)mxGetM(m_pm);
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return (index_t)mxGetN(m_pm);
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return (index_t)mxGetNumberOfElements(m_pm);
		}

		LMAT_ENSURE_INLINE const mwSize* dims() const
		{
			return mxGetDimensions(m_pm);
		}

	public:

		// Data access

		LMAT_ENSURE_INLINE double get_scalar() const
		{
			return mxGetScalar(m_pm);
		}

		template<typename T>
		LMAT_ENSURE_INLINE T scalar() const
		{
			return *((const T*)mxGetData(m_pm));
		}


		LMAT_ENSURE_INLINE const double *ptr_real() const
		{
			return mxGetPr(m_pm);
		}

		LMAT_ENSURE_INLINE const double *ptr_imag() const
		{
			return mxGetPi(m_pm);
		}

		template<typename T>
		LMAT_ENSURE_INLINE const T* data() const
		{
			return (const T*)mxGetData(m_pm);
		}

	public:

		// for objects, structs, and cells

		LMAT_ENSURE_INLINE
		const_marray get_property(index_t i, const char* name) const
		{
			return mxGetProperty(m_pm, (mwIndex)(i), name);
		}

		LMAT_ENSURE_INLINE
		const_marray get_field(index_t i, const char* name) const
		{
			return mxGetField(m_pm, (mwIndex)(i), name);
		}

		LMAT_ENSURE_INLINE
		int num_fields() const
		{
			return mxGetNumberOfFields(m_pm);
		}

		LMAT_ENSURE_INLINE
		const char* field_name_bynum(index_t fieldnum) const
		{
			return mxGetFieldNameByNumber(m_pm, fieldnum);
		}

		LMAT_ENSURE_INLINE
		int field_number(const char* fieldnum) const
		{
			return mxGetFieldNumber(m_pm, fieldnum);
		}

		LMAT_ENSURE_INLINE
		const_marray get_field_bynum(index_t i, int fieldnum) const
		{
			return mxGetFieldByNumber(m_pm, (mwIndex)(i), fieldnum);
		}

		LMAT_ENSURE_INLINE
		const_marray get_cell(index_t i) const
		{
			return mxGetCell(m_pm, (mwIndex)(i));
		}

	public:

		// for sparse matrices

		LMAT_ENSURE_INLINE index_t nzmax() const
		{
			return (index_t)mxGetNzmax(m_pm);
		}

		LMAT_ENSURE_INLINE const mwIndex* Ir() const
		{
			return mxGetIr(m_pm);
		}

		LMAT_ENSURE_INLINE const mwIndex* Jc() const
		{
			return mxGetJc(m_pm);
		}

	protected:
		mxArray *m_pm;

	}; // end class const_marray



	class marray : public const_marray
	{
	public:
		LMAT_ENSURE_INLINE marray(mxArray *pm)
		: const_marray(pm)
		{
		}

		LMAT_ENSURE_INLINE const mxArray *get() const
		{
			return this->m_pm;
		}

		LMAT_ENSURE_INLINE mxArray *get()
		{
			return this->m_pm;
		}

		LMAT_ENSURE_INLINE operator mxArray*()
		{
			return this->m_pm;
		}

		LMAT_ENSURE_INLINE void destroy()
		{
			mxDestroyArray(this->m_pm);
		}

	public:
		// Data access

		LMAT_ENSURE_INLINE const double *ptr_real() const
		{
			return mxGetPr(this->m_pm);
		}

		LMAT_ENSURE_INLINE double *ptr_real()
		{
			return mxGetPr(this->m_pm);
		}

		LMAT_ENSURE_INLINE const double *ptr_imag() const
		{
			return mxGetPi(this->m_pm);
		}

		LMAT_ENSURE_INLINE double *ptr_imag()
		{
			return mxGetPi(this->m_pm);
		}

		template<typename T>
		LMAT_ENSURE_INLINE const T* data() const
		{
			return (const T*)mxGetData(this->m_pm);
		}

		template<typename T>
		LMAT_ENSURE_INLINE T* data()
		{
			return (T*)mxGetData(this->m_pm);
		}

	public:

		// for objects, structs, and cells

		LMAT_ENSURE_INLINE
		void set_property(index_t i, const char* name, const_marray v)
		{
			mxSetProperty(this->m_pm, mwIndex(i), name, v.get());
		}

		LMAT_ENSURE_INLINE
		void set_field(index_t i, const char* name, marray v)
		{
			mxSetField(this->m_pm, mwIndex(i), name, v.get());
		}

		LMAT_ENSURE_INLINE
		void set_field_bynum(index_t i, int fieldnum, marray v)
		{
			mxSetFieldByNumber(this->m_pm, mwIndex(i), fieldnum, v.get());
		}

		LMAT_ENSURE_INLINE
		void set_cell(index_t i, marray v)
		{
			mxSetCell(this->m_pm, (mwIndex)(i), v.get());
		}

	public:
		// marray creating functions

		template<typename T>
		LMAT_ENSURE_INLINE
		static marray numeric_matrix(index_t m, index_t n, mxComplexity cplx=mxREAL)
		{
			const mxClassID cid = type_to_mclassid<T>::value;
			return mxCreateNumericMatrix(mwSize(m), mwSize(n), cid, cplx);
		}

		LMAT_ENSURE_INLINE
		static marray double_matrix(index_t m, index_t n, mxComplexity cplx=mxREAL)
		{
			return mxCreateDoubleMatrix(mwSize(m), mwSize(n), cplx);
		}

		LMAT_ENSURE_INLINE
		static marray logical_matrix(index_t m, index_t n)
		{
			return mxCreateLogicalMatrix(mwSize(m), mwSize(n));
		}

		template<typename T>
		LMAT_ENSURE_INLINE
		static marray numeric_cube(index_t m, index_t n, index_t k,
				mxComplexity cplx=mxREAL)
		{
			const mxClassID cid = type_to_mclassid<T>::value;

			mwSize dims[3] = {mwSize(m), mwSize(n), mwSize(k)};
			return mxCreateNumericArray(3, dims, cid, cplx);
		}

		LMAT_ENSURE_INLINE
		static marray double_cube(index_t m, index_t n, index_t k,
				mxComplexity cplx=mxREAL)
		{
			return numeric_cube<double>(m, n, k, cplx);
		}

		LMAT_ENSURE_INLINE
		static marray logical_cube(index_t m, index_t n, index_t k)
		{
			mwSize dims[3] = {mwSize(m), mwSize(n), mwSize(k)};
			return mxCreateLogicalArray(3, dims);
		}

		template<typename T>
		LMAT_ENSURE_INLINE
		static marray from_scalar(T x)
		{
			marray a = numeric_matrix<T>(1, 1);
			*(a.data<T>()) = x;
			return a;
		}

		LMAT_ENSURE_INLINE
		static marray from_double_scalar(double x)
		{
			return mxCreateDoubleScalar(x);
		}

		LMAT_ENSURE_INLINE
		static marray from_logical_scalar(bool v)
		{
			return mxCreateLogicalScalar(v);
		}

		LMAT_ENSURE_INLINE
		static marray from_string(const char* str)
		{
			return mxCreateString(str);
		}

		LMAT_ENSURE_INLINE
		static marray struct_matrix(int nfields, const char **fieldnames,
				index_t m = 1, index_t n = 1)
		{
			return mxCreateStructMatrix(mwSize(m), mwSize(n), nfields, fieldnames);
		}

		LMAT_ENSURE_INLINE
		static marray cell_matrix(index_t m, index_t n)
		{
			return mxCreateCellMatrix(mwSize(m), mwSize(n));
		}


	}; // end class marray


	LMAT_ENSURE_INLINE
	inline marray duplicate(const_marray a)
	{
		return mxDuplicateArray(a.get());
	}



} }

#endif /* MARRAY_H_ */
