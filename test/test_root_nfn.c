/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/root.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

void f(const double *in, double *out)
{
	double x = in[0], y = in[1], z = in[2];

	out[0] = x * x - y * y + z * z;
	out[1] = x + y - z;
	out[2] = x + y + z;
}

int main(void)
{
	double x[] = {1, 1, 1}; /* initial guess */
	double fx[3];
	double buffer[3 * (2 * 3 + 3)]; /* buffer length = n(2n + 3) */
	int ipvt[3];
	int code;
	int state = impf_root_nfn_newton(f, x, 3, buffer, ipvt, &code, 1e-8, 1e-8, 1000);

	printf("error = %i\n", code);
	assert(state == impf_EXIT_SUCCESS);
	assert(maxabs_arrd(fx, 3, 1) < 1e-8);
	return 0;
}
