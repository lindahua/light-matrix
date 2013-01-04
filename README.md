# Welcome to LightMatrix

**LightMatrix** is a C++ template library for matrix computation, which provides a unique combination of *design-time productivity* and *run-time performance*.
With this library, user can write codes using a set of friendly APIs (just like in MATLAB). These codes will then be transparently translated to highly optimized routines (at compile-time), through a meta-programming engine. 

## Overview

Relying a carefully designed core framework, LightMatrix has a series of nice properties. Here is a brief summary:

 - [LightMatrix is easy to use](#easy_to_use): one can perform computation easily using natural expressions
 - *LightMatrix is efficient:* high level codes will be mapped to tight loops of highly optimized SIMD instructions through meta-programming. No temporary matrices will be created in most cases.
 - *LightMatrix is convenient:* it supports various ways to access elements and sub parts of a matrix.
 - *LightMatrix is versatile:* it offers a broad range of functions out of box, which range from elementary functions, matrix reduction, to linear algebra.
 - *LightMatrix works nicely with other libraries:* it can directly work on external memories (without copying). Interfaces to work with STL vectors and MATLAB arrays are provided.
 - *LightMatrix supports kernel-based computation:* it has elegant syntax to define kernels and apply them to matrices in an element-wise manner. This provides a very efficient way for multi-input multi-output computation.
 - *LightMatrix is extensible:* one can easily introduce new matrix expressions and evaluation methods.
 - *LightMatrix is light-weight:* it is a pure header library, and it does not rely on third-party libraries such as BOOST.
 - *LightMatrix is reliable:* quality assured through extensive unit testing.
	
**Below are some specific examples/explanations to illustrate these features.**

<a id="easy_to_use">
### *LightMatrix* is easy to use
</a>

You can easily make matrices and perform computation using matrix expressions.

```c++	
using namespace lmat;

// generate a 2-by-3 matrix by inserting elements in row-major order
dense_matrix<double> a(2, 4, rm_({1., 2., 3., 4., 5., 6., 7., 8.}));
	
// generate a matrix by copying elements from a memory source 
double src[] = {1., 2., 3., 4., 5., 6., 7., 8.};
dense_matrix<double> b(2, 4, copy_from(src));
	
// do some computation
a += exp(b) * log(a + 1.0);
```