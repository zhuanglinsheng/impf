/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>
#include <math.h>
#include <string.h>

/* f(x, n; a, b, c) = (x1 + x2 + x3)^2 */
double f(const double *x, const int n)
{
	if (n < 3)
		return 0;
	else {
		double tmp = x[0] + x[1] + x[2];

		return tmp * tmp;
	}
}
/* f(x, n; a, b, c) = sin(x1) + sin(x2) + ... + sin(xn) */
double g(const double *x, const int n)
{
	double s = 0;
	int i;

	for (i = 0; i < n; i++)
		s += sin(x[i]);
	return s;
}
/* f(x, n; a, b, c) = cos(x1) + cos(x2) + ... + cos(xn) */
double h(const double *x, const int n)
{
	double s = 0;
	int i;

	for (i = 0; i < n; i++)
		s += cos(x[i]);
	return s;
}
/*           | (    x1  +     x2  + ... +     xn)^2 |
 * f(x, n) = |  sin(x1) + sin(x2) + ... + sin(xn)   |
 *           |  cos(x1) + cos(x2) + ... + cos(xn)   |
 */
void f2(const double *x, double *out)
{
	int i;

	memset(out, 0, 3 * sizeof(double));

	for (i = 0; i < 3; i++) {
		out[0] += x[i];
		out[1] += sin(x[i]);
		out[2] += cos(x[i]);
	}
	out[0] *= out[0];
}

int main(void)
{
	double x = 1. / 3.;
	double y = 2. / 3.;
	double z = 1.;
	double point[3];

	/* analytical solution */
	double fx = 2. * (x + y + z);
	double fy = 2. * (x + y + z);
	double fz = 2. * (x + y + z);
	double gx = cos(x);
	double gy = cos(y);
	double gz = cos(z);
	double hx = -sin(x);
	double hy = -sin(y);
	double hz = -sin(z);

	double Jac_n[9];
	double buffer[5 * 3]; /* buffer size = number of variables (n) */

	point[0] = x, point[1] = y, point[2] = z;

	impf_diff_1fn(f, point, 3, Jac_n, buffer);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);

	impf_diff_2fn(f, g, point, 3, Jac_n, buffer);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - gz) < 1e-8);

	impf_diff_3fn(f, g, h, point, 3, Jac_n, buffer);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - gz) < 1e-8);
	assert(__impf_ABS__(Jac_n[6] - hx) < 1e-6);
	assert(__impf_ABS__(Jac_n[7] - hy) < 1e-8);
	assert(__impf_ABS__(Jac_n[8] - hz) < 1e-8);

	impf_diff_mfn(f2, point, 3, Jac_n, 3, buffer);
	assert(__impf_ABS__(Jac_n[0] - fx) < 1e-8);
	assert(__impf_ABS__(Jac_n[1] - fy) < 1e-8);
	assert(__impf_ABS__(Jac_n[2] - fz) < 1e-8);
	assert(__impf_ABS__(Jac_n[3] - gx) < 1e-8);
	assert(__impf_ABS__(Jac_n[4] - gy) < 1e-8);
	assert(__impf_ABS__(Jac_n[5] - gz) < 1e-8);
	assert(__impf_ABS__(Jac_n[6] - hx) < 1e-6);
	assert(__impf_ABS__(Jac_n[7] - hy) < 1e-8);
	assert(__impf_ABS__(Jac_n[8] - hz) < 1e-8);
	return 0;
}
