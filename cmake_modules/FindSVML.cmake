# A CMake module to locate Intel SVML

# This script looks for two places:
#	- the environment variable ICC_LIBPATH
#	- the directory /opt/intel/lib

# Stage 1: find library directory

set(SVML_PATH $ENV{ICC_LIBPATH})

if (NOT SVML_PATH)
	# try to find at /opt/intel/lib
		
	if (EXISTS "/opt/intel/lib")
		set(SVML_PATH "/opt/intel/lib")
	endif (EXISTS "/opt/intel/lib")
endif (NOT SVML_PATH)


# Stage 2: find library files

if (SVML_PATH)

    if (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    	find_library(SVML svml_disp
    		HINTS ${SVML_PATH})
    	find_library(LIBIRC libirc 
    		HINTS ${SVML_PATH})
    		
    	set(SVML_LIBRARY
    	    ${SVML}
    	    ${LIBIRC})
    		
    else (${CMAKE_SYSTEM_NAME} MATCHES "Windows")
    	find_library(SVML svml 
    		HINTS ${SVML_PATH})
    		
    	set(SVML_LIBRARY
    	    ${SVML})
    	    
    endif (${CMAKE_SYSTEM_NAME} MATCHES "Windows")

endif (SVML_PATH)

# deal with QUIET and REQUIRED argument

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(SVML DEFAULT_MSG 
    SVML)
    
mark_as_advanced(SVML)
