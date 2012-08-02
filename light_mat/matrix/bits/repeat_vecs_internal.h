/**
 * @file repeat_vecs_internal.h
 *
 * Internal implementation for repeat_vecs
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef REPEAT_VECS_INTERNAL_H_
#define REPEAT_VECS_INTERNAL_H_

#include <light_mat/matrix/matrix_classes.h>

namespace lmat { namespace detail {

	template<class Col>
	struct repcol_cached_ewrapper
	{
		typedef typename matrix_traits<Col>::value_type T;
		typedef dense_matrix<T, ct_rows<Col>::value, 1> col_t;
		col_t m_col;

		LMAT_ENSURE_INLINE
		repcol_cached_ewrapper(const Col& col)
		: m_col(col) { }

		LMAT_ENSURE_INLINE
		index_t nrows() const { return m_col.nrows(); }

		LMAT_ENSURE_INLINE
		const T* data() const { return m_col.ptr_data(); }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t i) const { return m_col[i]; }

		LMAT_ENSURE_INLINE
		const col_t& ref() const { return m_col; }
	};

	template<class Col>
	struct repcol_ref_ewrapper
	{
		typedef typename matrix_traits<Col>::value_type T;
		typedef Col col_t;
		const col_t& m_col;

		LMAT_ENSURE_INLINE
		repcol_ref_ewrapper(const Col& col)
		: m_col(col) { }

		LMAT_ENSURE_INLINE
		index_t nrows() const { return m_col.nrows(); }

		LMAT_ENSURE_INLINE
		const T* data() const { return m_col.ptr_data(); }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t i) const { return m_col[i]; }

		LMAT_ENSURE_INLINE
		const col_t& ref() const { return m_col; }
	};

	template<class Row>
	struct reprow_cached_ewrapper
	{
		typedef typename matrix_traits<Row>::value_type T;
		typedef dense_matrix<T, 1, ct_cols<Row>::value> row_t;
		row_t m_row;

		LMAT_ENSURE_INLINE
		reprow_cached_ewrapper(const Row& row)
		: m_row(row) { }

		LMAT_ENSURE_INLINE
		index_t ncolumns() const { return m_row.ncolumns(); }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t j) const { return m_row[j]; }
	};

	template<class Row>
	struct reprow_ref_ewrapper
	{
		typedef typename matrix_traits<Row>::value_type T;
		typedef Row row_t;
		const row_t& m_row;

		LMAT_ENSURE_INLINE
		reprow_ref_ewrapper(const Row& row)
		: m_row(row) { }

		LMAT_ENSURE_INLINE
		index_t ncolumns() const { return m_row.ncolumns(); }

		LMAT_ENSURE_INLINE
		T operator[] (const index_t j) const { return m_row(0, j); }
	};


	template<class Col>
	struct repcol_ewrapper_map
	{
		typedef typename
				if_<is_dense_mat<Col>,
					repcol_ref_ewrapper<Col>,
					repcol_cached_ewrapper<Col> >::type type;
	};

	template<class Row>
	struct reprow_ewrapper_map
	{
		typedef typename
				if_<is_dense_mat<Row>,
					reprow_ref_ewrapper<Row>,
					reprow_cached_ewrapper<Row> >::type type;
	};

} }


#endif /* REPEAT_VECS_INTERNAL_H_ */
