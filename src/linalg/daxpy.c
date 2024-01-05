/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/linalg.h>

void impf_linalg_daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy)
{
	int i = 0, j = 0;

	assert(x != NULL);
	assert(y != NULL);

	while (i < n && j < n) {
		y[i] += a * x[j];
		i += incy;
		j += incx;
	}
}

void impf_linalg_ldaxpy(const int n, const long double a, const long double *x, const int incx,
			long double *y, const int incy)
{
	int i = 0, j = 0;

	assert(x != NULL);
	assert(y != NULL);

	while (i < n && j < n) {
		y[i] += a * x[j];
		i += incy;
		j += incx;
	}
}
