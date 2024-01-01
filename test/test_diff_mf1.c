/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>

double a = 1.0;
double b = 2.0;
double c = -3.0;

double f(const double x)
{
	return a * x * x + b * x + c;
}
double g(const double x)
{
	return a * x * x * x + b * x * x + c * x;
}
double h(const double x)
{
	return a * (x - b) * (x - c);
}

int main(void)
{
	double x = 1. / 3.;
	double Jac_n[3];

	/* analytical solution */
	double fx = 2. * a * x + b;
	double gx = 3. * a * x * x + 2. * b * x + c;
	double hx = a * (2. * x - b - c);

	impf_diff_1f1(f, x, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);

	impf_diff_2f1(f, g, x, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - gx) < 1e-8);

	impf_diff_3f1(f, g, h, x, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - hx) < 1e-8);
	return 0;
}
