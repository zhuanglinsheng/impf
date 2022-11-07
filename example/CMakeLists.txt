include_directories(${PROJECT_SOURCE_DIR}/include)
find_package(BLAS REQUIRED)

###############################################################################
# Drivers
###############################################################################

add_subdirectory(dts)
add_subdirectory(linalg)
add_subdirectory(ndiff)
add_subdirectory(root)

###############################################################################
# Examples
###############################################################################

add_executable(eg__heat2D eg__heat2D.c)
target_link_libraries(eg__heat2D impf m blas)
add_test(
    NAME eg__heat2D
    COMMAND eg__heat2D
)
