include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(test_linalg_dgels test_linalg_dgels.c)
target_link_libraries(test_linalg_dgels impf)
add_test(
	NAME test_linalg_dgels
	COMMAND test_linalg_dgels
)

add_executable(test_diff_mf1 test_diff_mf1.c)
target_link_libraries(test_diff_mf1 impf)
add_test(
	NAME test_diff_mf1
	COMMAND test_diff_mf1
)

add_executable(test_diff_mf2 test_diff_mf2.c)
target_link_libraries(test_diff_mf2 impf)
add_test(
	NAME test_diff_mf2
	COMMAND test_diff_mf2
)

add_executable(test_diff_mf3 test_diff_mf3.c)
target_link_libraries(test_diff_mf3 impf)
add_test(
	NAME test_diff_mf3
	COMMAND test_diff_mf3
)

add_executable(test_diff_mfn test_diff_mfn.c)
target_link_libraries(test_diff_mfn impf)
add_test(
	NAME test_diff_mfn
	COMMAND test_diff_mfn
)

add_executable(test_root_1f1 test_root_1f1.c)
target_link_libraries(test_root_1f1 impf m)
add_test(
	NAME test_root_1f1
	COMMAND test_root_1f1
)

add_executable(test_root_2f2 test_root_2f2.c)
target_link_libraries(test_root_2f2 impf m)
add_test(
	NAME test_root_2f2
	COMMAND test_root_2f2
)

add_executable(test_root_3f3 test_root_3f3.c)
target_link_libraries(test_root_3f3 impf m)
add_test(
	NAME test_root_3f3
	COMMAND test_root_3f3
)

add_executable(test_root_nfn test_root_nfn.c)
target_link_libraries(test_root_nfn impf m)
add_test(
	NAME test_root_nfn
	COMMAND test_root_nfn
)

add_executable(test_lp_simplex_1 test_lp_simplex_1.c)
target_link_libraries(test_lp_simplex_1 impf)
add_test(
	NAME test_lp_simplex_1
	COMMAND test_lp_simplex_1
)

add_executable(test_lp_simplex_2 test_lp_simplex_2.c)
target_link_libraries(test_lp_simplex_2 impf)
add_test(
	NAME test_lp_simplex_2
	COMMAND test_lp_simplex_2
)

add_executable(test_lp_simplex_3 test_lp_simplex_3.c)
target_link_libraries(test_lp_simplex_3 impf)
add_test(
	NAME test_lp_simplex_3
	COMMAND test_lp_simplex_3
)

add_executable(test_lp_simplex_4 test_lp_simplex_4.c)
target_link_libraries(test_lp_simplex_4 impf)
add_test(
	NAME test_lp_simplex_4
	COMMAND test_lp_simplex_4
)

add_executable(test_lp_simplex_5 test_lp_simplex_5.c)
target_link_libraries(test_lp_simplex_5 impf)
add_test(
	NAME test_lp_simplex_5
	COMMAND test_lp_simplex_5
)

add_executable(test_lp_simplex_6 test_lp_simplex_6.c)
target_link_libraries(test_lp_simplex_6 impf)
add_test(
	NAME test_lp_simplex_6
	COMMAND test_lp_simplex_6
)

add_executable(test_lp_simplex_netlib test_lp_simplex_netlib.c)
target_link_libraries(test_lp_simplex_netlib impf)
add_test(
	NAME test_lp_simplex_netlib
	COMMAND test_lp_simplex_netlib
)
