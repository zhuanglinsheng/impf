/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/root.h>
#include <math.h>
#include <stdio.h>

double a = 1.0;
double b = 2.0;
double c = -3.0;

double f(const double x)
{
	return a * x * x + b * x + c;
}

int main(void)
{
	double x, y;
	enum impf_ErrorCode code;
	int state;

	x = -b / (2. * a);
	y = x + 10.;
	printf("error = %i\n", code);
	state = impf_root_1f1_bisection(f, &x, &y, &code, 1e-8, 1e-8, 1000);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(f(x)) < 1e-8);

	x = 0.;
	state = impf_root_1f1_newton(f, &x, &code, 1e-8, 1e-8, 1000);
	printf("error = %i\n", code);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(f(x)) < 1e-8);

	x = -b / (2. * a);
	y = x + 10.;
	state = impf_root_1f1_secant(f, &x, &y, &code, 1e-8, 1e-8, 1000);
	printf("error = %i\n", code);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(f(x)) < 1e-8);

	return 0;
}
