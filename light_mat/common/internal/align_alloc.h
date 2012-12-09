/**
 * @file align_alloc.h
 *
 * @brief Implementation of aligned allocation
 *
 * @author Dahua Lin
 */

#ifdef _MSC_VER
#pragma once
#endif

#ifndef LIGHTMAT_ALIGN_ALLOC_H_
#define LIGHTMAT_ALIGN_ALLOC_H_

#include <light_mat/common/basic_defs.h>

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX
#include <stdlib.h>
#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32
#include <malloc.h>
#endif

namespace lmat { namespace internal {

#if LIGHTMAT_PLATFORM == LIGHTMAT_POSIX

	LMAT_ENSURE_INLINE
	inline void* aligned_allocate(size_t nbytes, unsigned int alignment)
	{
		char* p = 0;
		if (::posix_memalign((void**)(&p), alignment, nbytes) != 0)
		{
			throw std::bad_alloc();
		}
		return p;
	}

	LMAT_ENSURE_INLINE
	inline void aligned_release(void *p)
	{
		::free(p);
	}

#elif LIGHTMAT_PLATFORM == LIGHTMAT_WIN32

	LMAT_ENSURE_INLINE
	inline void* aligned_allocate(size_t nbytes, unsigned int alignment)
	{
		void* p = ::_aligned_malloc(nbytes, alignment));
		if (!p)
		{
			throw std::bad_alloc();
		}
		return p;
	}

	LMAT_ENSURE_INLINE
	inline void aligned_release(void *p)
	{
		::_aligned_free(p);
	}

#endif


} }

#endif /* ALIGN_ALLOC_H_ */
