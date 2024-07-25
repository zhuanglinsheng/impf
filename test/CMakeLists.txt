include_directories(${PROJECT_SOURCE_DIR}/include)

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

add_executable(test_lp_simplex_7 test_lp_simplex_7.c)
target_link_libraries(test_lp_simplex_7 impf)
add_test(
	NAME test_lp_simplex_7
	COMMAND test_lp_simplex_7
)

add_executable(test_lp_simplex_netlib test_lp_simplex_netlib.c)
target_link_libraries(test_lp_simplex_netlib impf)
add_test(
	NAME test_lp_simplex_netlib
	COMMAND test_lp_simplex_netlib
)
