add_executable(eg_root_impf_1f1_newton eg_root_impf_1f1_newton.c)
target_link_libraries(eg_root_impf_1f1_newton impf m blas)
add_test(
    NAME eg_root_impf_1f1_newton
    COMMAND eg_root_impf_1f1_newton
)

add_executable(eg_root_impf_2f2_newton eg_root_impf_2f2_newton.c)
target_link_libraries(eg_root_impf_2f2_newton impf m blas)
add_test(
    NAME eg_root_impf_2f2_newton
    COMMAND eg_root_impf_2f2_newton
)

add_executable(eg_root_impf_3f3_newton eg_root_impf_3f3_newton.c)
target_link_libraries(eg_root_impf_3f3_newton impf m blas)
add_test(
    NAME eg_root_impf_3f3_newton
    COMMAND eg_root_impf_3f3_newton
)
