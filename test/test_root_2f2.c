/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/root.h>
#include <math.h>
#include <stdio.h>

double f(const double x, const double y)
{
	return x * y + x * x - y;
}
double g(const double x, const double y)
{
	return x * x + y * y + x - y;
}

int main(void)
{
	/* numerical: initial guess */
	double x = 1.5;
	double y = 0.5;

	enum impf_ErrorCode code;
	int state = impf_root_2f2_newton(f, g, &x, &y, &code, 1e-8, 1e-8, 1000);

	printf("error = %i\n", code);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(f(x, y)) < 1e-8);
	assert(__impf_ABS__(g(x, y)) < 1e-8);
	return 0;
}
