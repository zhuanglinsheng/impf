add_executable(eg_root_newton_1f1 eg_root_newton_1f1.c)
target_link_libraries(eg_root_newton_1f1 impf m blas)
add_test(
    NAME eg_root_newton_1f1
    COMMAND eg_root_newton_1f1
)

add_executable(eg_root_newton_2f2 eg_root_newton_2f2.c)
target_link_libraries(eg_root_newton_2f2 impf m blas)
add_test(
    NAME eg_root_newton_2f2
    COMMAND eg_root_newton_2f2
)

add_executable(eg_root_newton_3f3 eg_root_newton_3f3.c)
target_link_libraries(eg_root_newton_3f3 impf m blas)
add_test(
    NAME eg_root_newton_3f3
    COMMAND eg_root_newton_3f3
)
