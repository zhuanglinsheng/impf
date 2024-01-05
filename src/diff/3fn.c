/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>

void impf_diff_3fn(double (*f)(const double*, const int),
		   double (*g)(const double*, const int),
		   double (*h)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer)
{
	impf_diff_1fn(f, x, n, out, buffer);
	impf_diff_1fn(g, x, n, out + n, buffer);
	impf_diff_1fn(h, x, n, out + n + n, buffer);
}
