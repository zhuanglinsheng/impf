/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>

void impf_diff_1f1(double (*f)(const double), const double x, double *out)
{
	double h, two_h, twelf_h;
	double p1, p2, p3, p4;

	assert(f != NULL);
	assert(out != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	p1 = x - two_h;
	p2 = x - h;
	p3 = x + h;
	p4 = x + two_h;
	*out = (f(p1) - f(p4) + 8. * (f(p3) - f(p2))) / twelf_h;
}

void impf_diff_1f2(double (*f)(const double, const double), const double x, const double y, double *out)
{
	double h, two_h, twelf_h;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;

	assert(f != NULL);
	assert(out != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	x1 = x - two_h;
	y1 = y - two_h;
	x2 = x - h;
	y2 = y - h;
	x3 = x + h;
	y3 = y + h;
	x4 = x + two_h;
	y4 = y + two_h;

	*out       = (f(x1, y) - f(x4, y) + 8. * (f(x3, y) - f(x2, y))) / twelf_h;
	*(out + 1) = (f(x, y1) - f(x, y4) + 8. * (f(x, y3) - f(x, y2))) / twelf_h;
}

void impf_diff_1f3(double (*f)(const double, const double, const double),
		   const double x, const double y, const double z, double *out)
{
	double h, two_h, twelf_h;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	double z1, z2, z3, z4;

	assert(f != NULL);
	assert(out != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	x1 = x - two_h;
	y1 = y - two_h;
	z1 = z - two_h;
	x2 = x - h;
	y2 = y - h;
	z2 = z - h;
	x3 = x + h;
	y3 = y + h;
	z3 = z + h;
	x4 = x + two_h;
	y4 = y + two_h;
	z4 = z + two_h;

	*out       = (f(x1, y, z) - f(x4, y, z) + 8. * (f(x3, y, z) - f(x2, y, z))) / twelf_h;
	*(out + 1) = (f(x, y1, z) - f(x, y4, z) + 8. * (f(x, y3, z) - f(x, y2, z))) / twelf_h;
	*(out + 2) = (f(x, y, z1) - f(x, y, z4) + 8. * (f(x, y, z3) - f(x, y, z2))) / twelf_h;
}

void impf_diff_1fn(double (*f)(const double*, const int), const double *x, const int n, double *out, double *buffer)
{
	double h, two_h, twelf_h;
	double tmp1, tmp2, tmp3, tmp4;
	int i;

	assert(f != NULL);
	assert(x != NULL);
	assert(out != NULL);
	assert(buffer != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	impf_memcpy(buffer, x, sizeof(double) * n);

	for (i = 0; i < n; i++) {
		*(buffer + i) -= two_h;
		tmp1 = f(buffer, n);     /* f(x-2h) */
		*(buffer + i) += h;
		tmp2 = f(buffer, n);     /* f(x-h)  */
		*(buffer + i) += two_h;
		tmp3 = f(buffer, n);     /* f(x+h)  */
		*(buffer + i) += h;
		tmp4 = f(buffer, n);     /* f(x+2h) */
		*(buffer + i) -= two_h;  /* restore buffer */
		*(out + i) = (tmp1 - tmp4 + 8. * (tmp3 - tmp2)) / twelf_h;
	}
}

void impf_diff_2f1(double (*f)(const double), double (*g)(const double), const double x, double *out)
{
	impf_diff_1f1(f, x, out);
	impf_diff_1f1(g, x, out + 1);
}

void impf_diff_2f2(double (*f)(const double, const double),
		   double (*g)(const double, const double),
		   const double x, const double y, double *out)
{
	impf_diff_1f2(f, x, y, out);
	impf_diff_1f2(g, x, y, out + 2);
}

void impf_diff_2f3(double (*f)(const double, const double, const double),
		   double (*g)(const double, const double, const double),
		   const double x, const double y, const double z, double *out)
{
	impf_diff_1f3(f, x, y, z, out);
	impf_diff_1f3(g, x, y, z, out + 3);
}

void impf_diff_2fn(double (*f)(const double*, const int), double (*g)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer)
{
	impf_diff_1fn(f, x, n, out, buffer);
	impf_diff_1fn(g, x, n, out + n, buffer);
}

void impf_diff_3f1(double (*f)(const double), double (*g)(const double), double (*h)(const double),
		   const double x, double *out)
{
	impf_diff_1f1(f, x, out);
	impf_diff_1f1(g, x, out + 1);
	impf_diff_1f1(h, x, out + 2);
}

void impf_diff_3f2(double (*f)(const double, const double),
		   double (*g)(const double, const double),
		   double (*h)(const double, const double),
		   const double x, const double y, double *out)
{
	impf_diff_1f2(f, x, y, out);
	impf_diff_1f2(g, x, y, out + 2);
	impf_diff_1f2(h, x, y, out + 4);
}

void impf_diff_3f3(double (*f)(const double, const double, const double),
		   double (*g)(const double, const double, const double),
		   double (*h)(const double, const double, const double),
		   const double x, const double y, const double z, double *out)
{
	impf_diff_1f3(f, x, y, z, out);
	impf_diff_1f3(g, x, y, z, out + 3);
	impf_diff_1f3(h, x, y, z, out + 6);
}

void impf_diff_3fn(double (*f)(const double*, const int),
		   double (*g)(const double*, const int),
		   double (*h)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer)
{
	impf_diff_1fn(f, x, n, out, buffer);
	impf_diff_1fn(g, x, n, out + n, buffer);
	impf_diff_1fn(h, x, n, out + n + n, buffer);
}

void impf_diff_mfn(void (*f)(const double*, double*),
		   const double *x, const int n, double *out, const int m, double *buffer)
{
	double h, two_h, twelf_h;
	double *point, *ret1, *ret2, *ret3, *ret4;
	int i, j;

	assert(f != NULL);
	assert(x != NULL);
	assert(out != NULL);
	assert(buffer != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	point = buffer;
	ret1 = buffer + n;
	ret2 = ret1 + m;
	ret3 = ret2 + m;
	ret4 = ret3 + m;

	impf_memcpy(point, x, sizeof(double) * n);

	for (j = 0; j < n; j++) {  /* Column j */
		*(point + j) -= two_h;
		f(point, ret1);  /* f(x_j - 2h) */
		*(point + j) += h;
		f(point, ret2);  /* f(x_j - h)  */
		*(point + j) += two_h;
		f(point, ret3);  /* f(x_j + h)  */
		*(point + j) += h;
		f(point, ret4);  /* f(x_j + 2h) */
		*(point + j) -= two_h;

		for (i = 0; i < m; i++) { /* Row i */
			double tmp1 = ret1[i] - ret4[i];
			double tmp2 = (ret3[i] - ret2[i]) * 8.;

			out[j + i * n] = (tmp1 + tmp2) / twelf_h;
		}
	}
}
