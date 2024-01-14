/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/root.h>
#include <math.h>
#include <stdio.h>

double f(const double x, const double y, const double z)
{
	return x * y * z + x * y - y * z - x * z;
}
double g(const double x, const double y, const double z)
{
	return sin(x) + sin(y) + sin(z);
}
double h(const double x, const double y, const double z)
{
	return cos(x) + cos(y) - cos(z);
}

int main(void)
{
	/* numerical: initial guess */
	double x = 1.8;
	double y = 0.5;
	double z = 0.5;

	int code;
	int state = impf_root_3f3_newton(f, g, h, &x, &y, &z, &code, 1e-8, 1e-8, 1000);

	printf("error = %i\n", code);
	assert(state == impf_EXIT_SUCCESS);
	assert(__impf_ABS__(f(x, y, z)) < 1e-8);
	assert(__impf_ABS__(g(x, y, z)) < 1e-8);
	assert(__impf_ABS__(h(x, y, z)) < 1e-8);
	return 0;
}
