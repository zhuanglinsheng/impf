/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/linalg.h>
#include <stdio.h>

int main(void)
{
	double A[] = {
		3., 5., 6.,
		4., 1., 0.,
		0., 0., 7.
	};
	double b[] = {5., 8., 9.};
	double x1[3], x2[3];
	int code;
	int ipvt[3];
	double buffer[100];
	int state;

	state = impf_linalg_dgels_3( \
			A[0], A[1], A[2], A[3], A[4], A[5], A[6], A[7], A[8],
			b[0], b[1], b[2], x1, &code);

	printf("Error code = %u\n", code);
	assert(state == EXIT_SUCCESS);
	state = impf_linalg_dgels_n(3, 1, A, 3, b, 1, ipvt, x2, buffer, &code);
	assert(state == EXIT_SUCCESS);
	assert(maxabs_arrd_gap(x1, x2, 3, 1) < 1e-8);
	return 0;
}
