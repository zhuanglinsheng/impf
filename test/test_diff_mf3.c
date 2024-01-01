/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>
#include <math.h>

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
	double x = 1. / 3.;
	double y = 2. / 3.;
	double z = 1.;
	double Jac_n[9];

	/* analytical solution */
	double fx = y * z + y - z;
	double fy = x * z + x - z;
	double fz = x * y - y - x;
	double gx = cos(x);
	double gy = cos(y);
	double gz = cos(z);
	double hx = -sin(x);
	double hy = -sin(y);
	double hz = sin(z);

	impf_diff_1f3(f, x, y, z, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);

	impf_diff_2f3(f, g, x, y, z, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - gz) < 1e-8);

	impf_diff_3f3(f, g, h, x, y, z, Jac_n);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - gz) < 1e-8);
	assert(__impf_ABS__(Jac_n[6] - hx) < 1e-8);
	assert(__impf_ABS__(Jac_n[7] - hy) < 1e-8);
	assert(__impf_ABS__(Jac_n[8] - hz) < 1e-8);
	return 0;
}
