add_executable(eg_diff_impf_grd_f2 eg_diff_impf_grd_f2.c)
target_link_libraries(eg_diff_impf_grd_f2 impf m blas)
add_test(
    NAME eg_diff_impf_grd_f2
    COMMAND eg_diff_impf_grd_f2
)

add_executable(eg_diff_impf_jacob eg_diff_impf_jacob.c)
target_link_libraries(eg_diff_impf_jacob impf m blas)
add_test(
    NAME eg_diff_impf_jacob
    COMMAND eg_diff_impf_jacob
)