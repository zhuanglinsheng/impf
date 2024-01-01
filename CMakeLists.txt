cmake_minimum_required(VERSION 3.12)
project(impf LANGUAGES C)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_STANDARD 90)
set(CMAKE_C_EXTENSIONS OFF)
if(MSVC)
  # Force to always compile with W4
  if(CMAKE_C_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  else()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /W4 -fsanitize=address")
  endif()
else()
  # Update if necessary
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-long-long -pedantic -fsanitize=address")
endif()
##
## Installation setting
##
## If user would like to choose their own default installation location,
## set the value of `CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT` to be
## `true` by switch the comments:
##
set(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT false)
##
## Then set the installation location below by changing the default value
## of `CMAKE_INSTALL_PREFIX`:
##
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX "$HOME/.local" CACHE PATH "..." FORCE)
endif()
set(CMAKE_INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")
##
## Third-party libraries
## add_compile_definitions(WITH_BLAS)

###############################################################################
# shared library
###############################################################################
set(SOURCES
	src/utils.c
	src/linalg_dscal.c
	src/linalg_daxpy.c
	src/linalg_dgels_L.c
	src/linalg_dgels_n.c
	src/linalg_dgels_2.c
	src/linalg_dgels_3.c
	src/diff_1f1.c
	src/diff_1f2.c
	src/diff_1f3.c
	src/diff_1fn.c
	src/diff_2f1.c
	src/diff_2f2.c
	src/diff_2f3.c
	src/diff_2fn.c
	src/diff_3f1.c
	src/diff_3f2.c
	src/diff_3f3.c
	src/diff_3fn.c
	src/diff_mfn.c
	src/root_1f1_bisection.c
	src/root_1f1_secant.c
	src/root_1f1_newton.c
	src/root_2f2_newton.c
	src/root_3f3_newton.c
	src/root_nfn_newton.c
	src/lp_readmps.c
	src/lp_simplex_std.c
	src/lp_simplex_gen.c
	src/lp_simplex_wrp.c
)
add_library(impf SHARED
	${SOURCES}
)
target_include_directories(impf PRIVATE
	${PROJECT_SOURCE_DIR}/include
)
target_link_libraries(impf)

###############################################################################
# Installation
###############################################################################
install(TARGETS impf LIBRARY
	DESTINATION lib
)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/impf
	DESTINATION include
)

###############################################################################
# Build Mode
#
# Choose the build mode either "Release" or "Debug" by
# cmake -DCMAKE_BUILD_TYPE=Release/Debug ..
###############################################################################
##
## (Default) The "Release" mode is for installation.
if (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Release"))
	set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wall -O2")
##
## The "Debug" mode is for unit testings.
elseif (CMAKE_BUILD_TYPE AND (CMAKE_BUILD_TYPE STREQUAL "Debug"))
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")
	enable_testing()
	add_subdirectory(test)
endif ()