/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>
#include <math.h>

double f(const double x, const double y)
{
	return x * y + x * x - y;
}
double g(const double x, const double y)
{
	return x * x + y * y + x - y;
}
double h(const double x, const double y)
{
	return sin(x) + sin(y);
}

int main(void)
{
	double x = 1. / 3.;
	double y = 2. / 3.;
	double Jac_n[6];

	/* analytical solution */
	double fx = y + 2. * x;
	double fy = x - 1.;
	double gx = 2. * x + 1;
	double gy = 2. * y - 1;
	double hx = cos(x);
	double hy = cos(y);

	impf_diff_1f2(f, x, y, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);

	impf_diff_2f2(f, g, x, y, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gy) < 1e-8);

	impf_diff_3f2(f, g, h, x, y, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - hx) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - hy) < 1e-8);
	return 0;
}
