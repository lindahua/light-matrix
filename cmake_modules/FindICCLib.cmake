# A CMake module to find Intel ICC libraries (typically shipped with ICC)

# This script looks for two places:
#	- the environment variable ICC_LIBPATH
#	- the directory /opt/intel/lib

# Stage 1: find library directory

set(ICCLIB_PATH $ENV{ICC_LIBPATH})

if (NOT ICCLIB_PATH)
	# try to find at /opt/intel/lib
		
	if (EXISTS "/opt/intel/lib")
		set(ICCLIB_PATH "/opt/intel/lib")
	endif (EXISTS "/opt/intel/lib")
endif (NOT ICCLIB_PATH)


# Stage 2: find library files

if (ICCLIB_PATH)

    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        find_library(INTEL_IMF libimf HINTS ${ICCLIB_PATH})
    	find_library(SVML svml_disp HINTS ${ICCLIB_PATH})
    	find_library(LIBIRC libirc HINTS ${ICCLIB_PATH})
    		
    	set(SVML_LIBRARY
    	    ${SVML}
    	    ${LIBIRC})
    		
    else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
        find_library(INTEL_IMF imf HINTS ${ICCLIB_PATH})
    	find_library(SVML svml HINTS ${ICCLIB_PATH})
    	find_library(LIBIRC irc HINTS ${ICCLIB_PATH})
    	
    		
    	set(SVML_LIBRARY
    	    ${SVML})
    	    
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

endif (ICCLIB_PATH)

# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(ICCLib DEFAULT_MSG 
    INTEL_IMF SVML LIBIRC)
    
mark_as_advanced(ICCLib INTEL_IMF SVML LIBIRC)

