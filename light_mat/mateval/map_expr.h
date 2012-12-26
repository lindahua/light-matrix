/**
 * @file map_expr.h
 *
 * @brief Element-wise mapping expressions
 *
 * @author Dahua Lin
 */

#ifndef LIGHTMAT_MAP_EXPR_H_
#define LIGHTMAT_MAP_EXPR_H_

#include "internal/map_expr_internal.h"

namespace lmat
{

	// forward declarations

	template<typename... Args> class map_expr;


	/********************************************
	 *
	 *  Map expression classes
	 *
	 ********************************************/

	template<typename FTag, typename... Args>
	struct matrix_traits<map_expr<FTag, Args...> >
	{
		typedef typename internal::map_expr_helper<Args...>::type helper_t;

		static const int num_dimensions = 2;
		static const int ct_num_rows = helper_t::ct_nrows;
		static const int ct_num_cols = helper_t::ct_ncols;

		static const bool is_readonly = true;

		typedef typename helper_t::shape_type shape_type;

		typedef typename internal::map_expr_value<FTag, Args...>::type value_type;
		typedef typename helper_t::domain domain;
	};

	namespace meta
	{
		template<typename FTag, typename... Args>
		struct supports_ewise_access<map_expr<FTag, Args...> >
		{
			static const bool value = true;
		};
	}

	template<typename FTag, typename Arg1>
	class map_expr<FTag, Arg1>
	: public IMatrixXpr<map_expr<FTag, Arg1>,
	  typename internal::map_expr_value<FTag, Arg1>::type>
	{
		typedef typename internal::map_expr_helper<Arg1>::type helper_t;
		typedef typename helper_t::shape_type shape_type;

	public:
		LMAT_ENSURE_INLINE
		map_expr(FTag, const Arg1& a1)
		: m_arg1(a1) { }

		LMAT_ENSURE_INLINE FTag tag() const
		{
			return FTag();
		}

		LMAT_ENSURE_INLINE const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_arg1.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_arg1.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_arg1.nelems();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_arg1.shape();
		}

	private:
		const Arg1& m_arg1;
	};


	template<typename FTag, typename Arg1, typename Arg2>
	class map_expr<FTag, Arg1, Arg2>
	: public IMatrixXpr<map_expr<FTag, Arg1, Arg2>,
	  typename internal::map_expr_value<FTag, Arg1, Arg2>::type>
	{
		typedef typename internal::map_expr_helper<Arg1, Arg2>::type helper_t;
		typedef typename helper_t::shape_type shape_type;

	public:
		LMAT_ENSURE_INLINE
		map_expr(FTag, const Arg1& a1, const Arg2& a2)
		: m_shape(helper_t::get_shape(a1, a2))
		, m_arg1(a1), m_arg2(a2) { }

		LMAT_ENSURE_INLINE FTag tag() const
		{
			return FTag();
		}

		LMAT_ENSURE_INLINE const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE const Arg2& arg2() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
	};


	template<typename FTag, typename Arg1, typename Arg2, typename Arg3>
	class map_expr<FTag, Arg1, Arg2, Arg3>
	: public IMatrixXpr<map_expr<FTag, Arg1, Arg2, Arg3>,
	  typename internal::map_expr_value<FTag, Arg1, Arg2, Arg3>::type>
	{
		typedef typename internal::map_expr_helper<Arg1, Arg2, Arg3>::type helper_t;
		typedef typename helper_t::shape_type shape_type;

	public:
		LMAT_ENSURE_INLINE
		map_expr(FTag, const Arg1& a1, const Arg2& a2, const Arg3& a3)
		: m_shape(helper_t::get_shape(a1, a2, a3))
		, m_arg1(a1), m_arg2(a2), m_arg3(a3) { }

		LMAT_ENSURE_INLINE FTag tag() const
		{
			return FTag();
		}

		LMAT_ENSURE_INLINE const Arg1& arg1() const
		{
			return m_arg1;
		}

		LMAT_ENSURE_INLINE const Arg2& arg2() const
		{
			return m_arg2;
		}

		LMAT_ENSURE_INLINE const Arg3& arg3() const
		{
			return m_arg3;
		}

		LMAT_ENSURE_INLINE index_t nrows() const
		{
			return m_shape.nrows();
		}

		LMAT_ENSURE_INLINE index_t ncolumns() const
		{
			return m_shape.ncolumns();
		}

		LMAT_ENSURE_INLINE index_t nelems() const
		{
			return m_shape.nelems();
		}

		LMAT_ENSURE_INLINE shape_type shape() const
		{
			return m_shape;
		}

	private:
		shape_type m_shape;
		const Arg1& m_arg1;
		const Arg2& m_arg2;
		const Arg3& m_arg3;
	};


	/********************************************
	 *
	 *  Expression construction functions
	 *
	 ********************************************/

	template<typename FTag, class A1, typename T1>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A1>,
		map_expr<FTag, A1> >::type
	make_map_expr(const FTag& ftag, const IMatrixXpr<A1, T1>& a1)
	{
		return map_expr<FTag, A1>(ftag, a1.derived());
	}


	template<typename FTag, class A1, typename T1, typename T2>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A1>,
		map_expr<FTag, A1, T2> >::type
	make_map_expr_fix2(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, const T2& a2)
	{
		return map_expr<FTag, A1, T2>(ftag, a1.derived(), a2);
	}

	template<typename FTag, typename T1, class A2, typename T2>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A2>,
		map_expr<FTag, T1, A2> >::type
	make_map_expr_fix1(const FTag& ftag, const T1& a1, const IMatrixXpr<A2, T2>& a2)
	{
		return map_expr<FTag, T1, A2>(ftag, a1, a2.derived());
	}

	template<typename FTag, typename A1, typename T1, typename T2, class A2>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::and_<meta::supports_ewise_access<A1>, meta::supports_ewise_access<A2> >,
		map_expr<FTag, A1, A2> >::type
	make_map_expr(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, const IMatrixXpr<A2, T2>& a2)
	{
		return map_expr<FTag, A1, A2>(ftag, a1.derived(), a2.derived());
	}


	template<typename FTag, class A1, typename T1, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A1>,
		map_expr<FTag, A1, T2, T3> >::type
	make_map_expr_fix23(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, const T2& a2, const T3& a3)
	{
		return map_expr<FTag, A1, T2, T3>(ftag, a1.derived(), a2, a3);
	}

	template<typename FTag, typename T1, class A2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A2>,
		map_expr<FTag, T1, A2, T3> >::type
	make_map_expr_fix13(const FTag& ftag, T1& a1, const IMatrixXpr<A2, T2>& a2, const T3& a3)
	{
		return map_expr<FTag, T1, A2, T3>(ftag, a1, a2.derived(), a3);
	}

	template<typename FTag, typename T1, typename T2, class A3, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::supports_ewise_access<A3>,
		map_expr<FTag, T1, T2, A3> >::type
	make_map_expr_fix12(const FTag& ftag, const T1& a1, const T2& a2, const IMatrixXpr<A3, T3>& a3)
	{
		return map_expr<FTag, T1, T2, A3>(ftag, a1, a2, a3.derived());
	}

	template<typename FTag, class A1, typename T1, class A2, typename T2, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::and_<meta::supports_ewise_access<A1>, meta::supports_ewise_access<A2> >,
		map_expr<FTag, A1, A2, T3> >::type
	make_map_expr_fix3(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, const IMatrixXpr<A2, T2>& a2, const T3& a3)
	{
		return map_expr<FTag, A1, A2, T3>(ftag, a1.derived(), a2.derived(), a3);
	}

	template<typename FTag, class A1, typename T1, typename T2, class A3, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::and_<meta::supports_ewise_access<A1>, meta::supports_ewise_access<A3> >,
		map_expr<FTag, A1, T2, A3> >::type
	make_map_expr_fix2(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, T2& a2, const IMatrixXpr<A3, T3>& a3)
	{
		return map_expr<FTag, A1, T2, A3>(ftag, a1.derived(), a2, a3.derived());
	}

	template<typename FTag, typename T1, class A2, typename T2, class A3, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if<
		meta::and_<meta::supports_ewise_access<A2>, meta::supports_ewise_access<A3> >,
		map_expr<FTag, T1, A2, A3> >::type
	make_map_expr_fix1(const FTag& ftag, const T1& a1, const IMatrixXpr<A2, T2>& a2, const IMatrixXpr<A3, T3>& a3)
	{
		return map_expr<FTag, T1, A2, A3>(ftag, a1, a2.derived(), a3.derived());
	}

	template<typename FTag, class A1, typename T1, class A2, typename T2, class A3, typename T3>
	LMAT_ENSURE_INLINE
	inline typename meta::enable_if_c<
		meta::supports_ewise_access<A1>::value &&
		meta::supports_ewise_access<A2>::value &&
		meta::supports_ewise_access<A3>::value,
		map_expr<FTag, A1, A2, A3> >::type
	make_map_expr(const FTag& ftag, const IMatrixXpr<A1, T1>& a1, const IMatrixXpr<A2, T2>& a2, const IMatrixXpr<A3, T3>& a3)
	{
		return map_expr<FTag, A1, A2, A3>(ftag, a1.derived(), a2.derived(), a3.derived());
	}


	/********************************************
	 *
	 *  Accessor classes
	 *
	 ********************************************/

	namespace internal
	{
		template<typename FTag, typename Arg1, typename U>
		struct vec_reader_map<map_expr<FTag, Arg1>, U>
		{
			typedef map_expr<FTag, Arg1> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1>::type fun_type;

			typedef typename internal::arg_vec_reader_map<Arg1, U>::type arg1_rd_t;
			typedef map_vec_reader<fun_type, U, arg1_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_vec_reader_map<Arg1, U>::get(expr.arg1()) );
			}
		};

		template<typename FTag, typename Arg1, typename Arg2, typename U>
		struct vec_reader_map<map_expr<FTag, Arg1, Arg2>, U>
		{
			typedef map_expr<FTag, Arg1, Arg2> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1, Arg2>::type fun_type;

			typedef typename internal::arg_vec_reader_map<Arg1, U>::type arg1_rd_t;
			typedef typename internal::arg_vec_reader_map<Arg2, U>::type arg2_rd_t;
			typedef map_vec_reader<fun_type, U, arg1_rd_t, arg2_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_vec_reader_map<Arg1, U>::get(expr.arg1()),
						internal::arg_vec_reader_map<Arg2, U>::get(expr.arg2()) );
			}
		};

		template<typename FTag, typename Arg1, typename Arg2, typename Arg3, typename U>
		struct vec_reader_map<map_expr<FTag, Arg1, Arg2, Arg3>, U>
		{
			typedef map_expr<FTag, Arg1, Arg2, Arg3> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1, Arg2, Arg3>::type fun_type;

			typedef typename internal::arg_vec_reader_map<Arg1, U>::type arg1_rd_t;
			typedef typename internal::arg_vec_reader_map<Arg2, U>::type arg2_rd_t;
			typedef typename internal::arg_vec_reader_map<Arg3, U>::type arg3_rd_t;
			typedef map_vec_reader<fun_type, U, arg1_rd_t, arg2_rd_t, arg3_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_vec_reader_map<Arg1, U>::get(expr.arg1()),
						internal::arg_vec_reader_map<Arg2, U>::get(expr.arg2()),
						internal::arg_vec_reader_map<Arg3, U>::get(expr.arg3()) );
			}
		};


		template<typename FTag, typename Arg1, typename U>
		struct multicol_reader_map<map_expr<FTag, Arg1>, U>
		{
			typedef map_expr<FTag, Arg1> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1>::type fun_type;

			typedef typename internal::arg_multicol_reader_map<Arg1, U>::type arg1_rd_t;
			typedef map_multicol_reader<fun_type, U, arg1_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_multicol_reader_map<Arg1, U>::get(expr.arg1()) );
			}
		};

		template<typename FTag, typename Arg1, typename Arg2, typename U>
		struct multicol_reader_map<map_expr<FTag, Arg1, Arg2>, U>
		{
			typedef map_expr<FTag, Arg1, Arg2> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1, Arg2>::type fun_type;

			typedef typename internal::arg_multicol_reader_map<Arg1, U>::type arg1_rd_t;
			typedef typename internal::arg_multicol_reader_map<Arg2, U>::type arg2_rd_t;
			typedef map_multicol_reader<fun_type, U, arg1_rd_t, arg2_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_multicol_reader_map<Arg1, U>::get(expr.arg1()),
						internal::arg_multicol_reader_map<Arg2, U>::get(expr.arg2()) );
			}
		};

		template<typename FTag, typename Arg1, typename Arg2, typename Arg3, typename U>
		struct multicol_reader_map<map_expr<FTag, Arg1, Arg2, Arg3>, U>
		{
			typedef map_expr<FTag, Arg1, Arg2, Arg3> expr_type;
			typedef typename internal::map_expr_fun<FTag, U, Arg1, Arg2, Arg3>::type fun_type;

			typedef typename internal::arg_multicol_reader_map<Arg1, U>::type arg1_rd_t;
			typedef typename internal::arg_multicol_reader_map<Arg2, U>::type arg2_rd_t;
			typedef typename internal::arg_multicol_reader_map<Arg3, U>::type arg3_rd_t;
			typedef map_multicol_reader<fun_type, U, arg1_rd_t, arg2_rd_t, arg3_rd_t> type;

			LMAT_ENSURE_INLINE
			static type get(const expr_type& expr)
			{
				return type(fun_type(), U(),
						internal::arg_multicol_reader_map<Arg1, U>::get(expr.arg1()),
						internal::arg_multicol_reader_map<Arg2, U>::get(expr.arg2()),
						internal::arg_multicol_reader_map<Arg3, U>::get(expr.arg3()) );
			}
		};
	}


	/********************************************
	 *
	 *  Evaluation
	 *
	 ********************************************/

	namespace internal
	{
		template<typename FTag, typename Arg1>
		struct prefers_linear<map_expr<FTag, Arg1> >
		{
			static const bool value = arg_prefers_linear<Arg1>::value;
		};

		template<typename FTag, typename Arg1, typename Arg2>
		struct prefers_linear<map_expr<FTag, Arg1, Arg2> >
		{
			static const bool value =
					arg_prefers_linear<Arg1>::value &&
					arg_prefers_linear<Arg2>::value;
		};

		template<typename FTag, typename Arg1, typename Arg2, typename Arg3>
		struct prefers_linear<map_expr<FTag, Arg1, Arg2, Arg3> >
		{
			static const bool value =
					arg_prefers_linear<Arg1>::value &&
					arg_prefers_linear<Arg2>::value &&
					arg_prefers_linear<Arg3>::value;
		};


		template<typename FTag, typename Arg1,
			typename T, typename Kind, bool IsLinear>
		struct prefers_simd<map_expr<FTag, Arg1>, T, Kind, IsLinear>
		{
			static const bool value =
					has_simd_support<FTag, T, Kind>::value &&
					arg_prefers_simd<Arg1, Kind, IsLinear>::value;
		};

		template<typename FTag, typename Arg1, typename Arg2,
			typename T, typename Kind, bool IsLinear>
		struct prefers_simd<map_expr<FTag, Arg1, Arg2>, T, Kind, IsLinear>
		{
			static const bool value =
					has_simd_support<FTag, T, Kind>::value &&
					arg_prefers_simd<Arg1, Kind, IsLinear>::value &&
					arg_prefers_simd<Arg2, Kind, IsLinear>::value;
		};

		template<typename FTag, typename Arg1, typename Arg2, typename Arg3,
			typename T, typename Kind, bool IsLinear>
		struct prefers_simd<map_expr<FTag, Arg1, Arg2, Arg3>, T, Kind, IsLinear>
		{
			static const bool value =
					has_simd_support<FTag, T, Kind>::value &&
					arg_prefers_simd<Arg1, Kind, IsLinear>::value &&
					arg_prefers_simd<Arg2, Kind, IsLinear>::value &&
					arg_prefers_simd<Arg3, Kind, IsLinear>::value;
		};
	}


	template<typename U>
	struct linear_map_policy
	{
		template<typename T, class Expr, class DMat>
		LMAT_ENSURE_INLINE
		static void eval(const IMatrixXpr<Expr, T>& expr, IRegularMatrix<DMat, T>& dmat)
		{
			const Expr& s = expr.derived();
			DMat& d = dmat.derived();

			ewise(copy_kernel<T>(), U())(common_shape(s, d), in_(s), out_(d));
		}
	};


	template<typename U>
	struct percol_map_policy
	{
		template<typename T, class Expr, class DMat>
		LMAT_ENSURE_INLINE
		static void eval(const IMatrixXpr<Expr, T>& expr, IRegularMatrix<DMat, T>& dmat)
		{
			const Expr& s = expr.derived();
			DMat& d = dmat.derived();

			percol(ewise(copy_kernel<T>(), U()), common_shape(s, d), in_(s), out_(d));
		}
	};


	template<typename FTag, typename... Args, class DMat>
	LMAT_ENSURE_INLINE
	inline void evaluate(const map_expr<FTag, Args...>& sexpr,
			IRegularMatrix<DMat, typename internal::map_expr_value<FTag, Args...>::type>& dmat)
	{
		typedef map_expr<FTag, Args...> expr_type;
		typedef typename internal::preferred_map_policy<expr_type, DMat>::type policy_t;

		policy_t::eval(sexpr, dmat.derived());
	}


}

#endif /* EWISE_EXPR_H_ */
