add_executable(eg_linalg_lsolve_2A2_2b eg_linalg_lsolve_2A2_2b.c)
target_link_libraries(eg_linalg_lsolve_2A2_2b impf m blas)
add_test(
    NAME eg_linalg_lsolve_2A2_2b
    COMMAND eg_linalg_lsolve_2A2_2b
)

add_executable(eg_linalg_lsolve_3A3_3b eg_linalg_lsolve_3A3_3b.c)
target_link_libraries(eg_linalg_lsolve_3A3_3b impf m blas)
add_test(
    NAME eg_linalg_lsolve_3A3_3b
    COMMAND eg_linalg_lsolve_3A3_3b
)