cmake_minimum_required(VERSION 3.0)
project(impf LANGUAGES C)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)


###############################################################################
# shared library
###############################################################################
set(SOURCES
    src/impf_blas_datastructs.c
    src/impf_blas_distances.c
    src/impf_blas_operations.c
    src/impf_blas_multiplications.c
    src/impf_linalg.c
    src/impf_ndiff.c
    src/impf_root.c
    src/impf_approx_poly.c
    src/impf_approx_fem.c
)
add_library(impf SHARED
    ${SOURCES}
)
target_include_directories(impf PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

###############################################################################
# Link to BLAS (for non-Linux users)
###############################################################################
find_package(BLAS REQUIRED)
target_link_libraries(impf m blas)

###############################################################################
# Installation
###############################################################################
## If user would like to choose their own default installation location, 
## set the value of `CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT` to be 
## `true` by switch the comments:  
##
set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT false)
# set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT true)
##
## Then set the installation location below by changing the default value 
## of `CMAKE_INSTALL_PREFIX`:
##
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "~/.local" CACHE PATH "..." FORCE)
endif()
## 
## Installations: Shared library and headers
## 
install(TARGETS impf LIBRARY 
    DESTINATION lib
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/impf 
    DESTINATION include
)


###############################################################################
# Build Mode
#
# Choose the build mode either "Release", "Debug" or "Doc" by 
# cmake -DCMAKE_BUILD_TYPE=Release/Debug/Doc ..
###############################################################################
##
## (Default) The "Release" mode is for installation, .
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -O2")
##
## The "Doc" mode is for documentations generation. 
elseif (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Doc"))
    find_package(Doxygen)
    if (DOXYGEN_FOUND)
        add_custom_target(doc_doxygen ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${PROJECT_SOURCE_DIR}/Doxyfile
            WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
            COMMENT "Generating API documentation with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message("Doxygen need to be installed to generate the documentation")
    endif (DOXYGEN_FOUND)
##
## The "Debug" mode is for unit testings.
elseif (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
    enable_testing()
    add_subdirectory(example)
endif ()

