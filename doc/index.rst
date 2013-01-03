.. LightMatrix documentation master file

Welcome to LightMatrix
=========================

**LightMatrix** is a C++ template library for matrix computation, which provides a unique combination of *design-time productivity* and *run-time performance*.
With this library, user can write codes using a set of friendly APIs (just like in MATLAB). These codes will then be transparently translated to highly optimized routines (at compile-time), through a meta-programming engine. 

Overview
----------

* **LightMatrix** is easy to use.

	You can easily make matrices and perform computation using matrix expressions.

	.. code-block:: c++

		using namespace lmat;

		dense_matrix<double> a = rand<double>(3, 5);
		dense_matrix<double> b = rand<double>(3, 5);
		a += exp(b) * log(a + 1.0);

* **LightMatrix** is efficient.

	When AVX is enabled, the last statement above will be translated transparently at compile time to something like the following: (*note: it also supports SSE*)

	.. code-block:: c++

		v1 = _mm256_set1_pd(1.0);
		for (index_t i = 0; i < a.nelems(); i += 4) 
		{
			va = _mm256_load_pd(a.ptr_data() + i);
			vb = _mm256_load_pd(b.ptr_data() + i);
			t1 = __svml_exp4(vb);
			t2 = __svml_log4(_mm256_add_pd(va, v1));
			r = _mm256_add_pd(t1, t2);
			_mm256_store_pd(a.ptr_data() + i, r); 
		}

	The performance is comparable to hand-optimized SIMD codes. (Node: ``__svml_exp4`` and ``__svml_log4`` are routines to compute exp and log on AVX vectors.)
	When explicitly requested, it can also generate multi-threaded codes using OpenMP.




Documentation
===============

.. toctree::
   :maxdepth: 2

   get_started
   basics
   advanced
   dev_guide



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

