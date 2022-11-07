add_executable(eg_dts_impf_t_mat eg_dts_impf_t_mat.c)
target_link_libraries(eg_dts_impf_t_mat impf m blas)
add_test(
    NAME eg_dts_impf_t_mat
    COMMAND eg_dts_impf_t_mat
)

add_executable(eg_dts_impf_t_vec eg_dts_impf_t_vec.c)
target_link_libraries(eg_dts_impf_t_vec impf m blas)
add_test(
    NAME eg_dts_impf_t_vec
    COMMAND eg_dts_impf_t_vec
)

add_executable(eg_dts_impf_mat_zeros eg_dts_impf_mat_zeros.c)
target_link_libraries(eg_dts_impf_mat_zeros impf m blas)
add_test(
    NAME eg_dts_impf_mat_zeros
    COMMAND eg_dts_impf_mat_zeros
)

add_executable(eg_dts_impf_mat_ones eg_dts_impf_mat_ones.c)
target_link_libraries(eg_dts_impf_mat_ones impf m blas)
add_test(
    NAME eg_dts_impf_mat_ones
    COMMAND eg_dts_impf_mat_ones
)

add_executable(eg_dts_impf_mat_eyes eg_dts_impf_mat_eyes.c)
target_link_libraries(eg_dts_impf_mat_eyes impf m blas)
add_test(
    NAME eg_dts_impf_mat_eyes
    COMMAND eg_dts_impf_mat_eyes
)

add_executable(eg_dts_impf_mat_transpose eg_dts_impf_mat_transpose.c)
target_link_libraries(eg_dts_impf_mat_transpose impf m blas)
add_test(
    NAME eg_dts_impf_mat_transpose
    COMMAND eg_dts_impf_mat_transpose
)

add_executable(eg_dts_impf_mat_transmajor eg_dts_impf_mat_transmajor.c)
target_link_libraries(eg_dts_impf_mat_transmajor impf m blas)
add_test(
    NAME eg_dts_impf_mat_transmajor
    COMMAND eg_dts_impf_mat_transmajor
)

add_executable(eg_dts_impf_mat_transmajor_inplace eg_dts_impf_mat_transmajor_inplace.c)
target_link_libraries(eg_dts_impf_mat_transmajor_inplace impf m blas)
add_test(
    NAME eg_dts_impf_mat_transmajor_inplace
    COMMAND eg_dts_impf_mat_transmajor_inplace
)

add_executable(eg_dts_impf_submat_subtract eg_dts_impf_submat_subtract.c)
target_link_libraries(eg_dts_impf_submat_subtract impf m blas)
add_test(
    NAME eg_dts_impf_submat_subtract
    COMMAND eg_dts_impf_submat_subtract
)

add_executable(eg_dts_impf_df_norms eg_dts_impf_df_norms.c)
target_link_libraries(eg_dts_impf_df_norms impf m blas)
add_test(
    NAME eg_dts_impf_df_norms
    COMMAND eg_dts_impf_df_norms
)

add_executable(eg_dts_impf_df_distances eg_dts_impf_df_distances.c)
target_link_libraries(eg_dts_impf_df_distances impf m blas)
add_test(
    NAME eg_dts_impf_df_distances
    COMMAND eg_dts_impf_df_distances
)