# a simple cmake script to locate Intel Math Kernel Library

# This script looks for two places:
#	- the environment variable MKLROOT
#	- the directory /opt/intel/mkl


# Stage 1: find the root directory

set(MKLROOT_PATH $ENV{MKLROOT})

if (NOT MKLROOT_PATH)
	# try to find at /opt/intel/mkl
		
	if (EXISTS "/opt/intel/mkl")
		set(MKLROOT_PATH "/opt/intel/mkl")
	endif (EXISTS "/opt/intel/mkl")
endif (NOT MKLROOT_PATH)


# Stage 2: find include path and libraries
	
if (MKLROOT_PATH)
	# root-path found
	
	set(EXPECT_MKL_INCPATH "${MKLROOT_PATH}/include")
	set(EXPECT_MKL_LIBPATH "${MKLROOT_PATH}/lib")
	
	if (IS_DIRECTORY ${EXPECT_MKL_INCPATH})
		set(MKL_INCLUDE_DIR ${EXPECT_MKL_INCPATH})
	endif (IS_DIRECTORY ${EXPECT_MKL_INCPATH})
	
	if (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
		set(MKL_LIBRARY ${EXPECT_MKL_LIBPATH})
	endif (IS_DIRECTORY ${EXPECT_MKL_LIBPATH})
	
endif (MKLROOT_PATH)

# Stage 3: linking libraries

if (MKL_LIBRARY)

	if (CMAKE_SYSTEM_NAME MATCHES "Darwin")
		if (CMAKE_SIZEOF_VOID_P MATCHES 8)
			set(MKL_LINK_FLAGS "-L${MKL_LIBRARY} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm")
		else (CMAKE_SIZEOF_VOID_P MATCHES 8)
			set(MKL_LINK_FLAGS "-L${MKL_LIBRARY} -lmkl_intel -lmkl_intel_thread -lmkl_core -openmp -lpthread -lm")
		endif (CMAKE_SIZEOF_VOID_P MATCHES 8)
	endif (CMAKE_SYSTEM_NAME MATCHES "Darwin")
	
	if (CMAKE_SYSTEM_NAME MATCHES "Linux")
		if (CMAKE_SIZEOF_VOID_P MATCHES 8)
			set(MKL_LINK_FLAGS "-L${MKL_LIBRARY} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
		else (CMAKE_SIZEOF_VOID_P MATCHES 8)
			set(MKL_LINK_FLAGS "-L${MKL_LIBRARY} -lmkl_intel -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm")
		endif (CMAKE_SIZEOF_VOID_P MATCHES 8)
	endif (CMAKE_SYSTEM_NAME MATCHES "Linux")

endif (MKL_LIBRARY)
	
# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(MKL DEFAULT_MSG MKL_LIBRARY MKL_INCLUDE_DIR)
mark_as_advanced(MKL_LIBRARY MKL_INCLUDE_DIR)


