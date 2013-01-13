# The module to locate AMD LibM

# This script looks for:
#	- the environment variable AMDLIBM_PATH
#

# Stage 1: find library directory

set(LIBM_PATH "$ENV{AMDLIBM_PATH}")

# Stage 2: find library file

if (LIBM_PATH)

    set(LIBM_LIBPATH "$ENV{AMDLIBM_PATH}/lib/dynamic")

    find_library(LIBM amdlibm
        HINTS ${LIBM_LIBPATH})

    set(LIBM_LIBRARY ${LIBM})

endif (LIBM_PATH)

# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(LIBM DEFAULT_MSG 
    LIBM)
    
mark_as_advanced(LIBM)


