/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>

void impf_diff_3f2(double (*f)(const double, const double),
		   double (*g)(const double, const double),
		   double (*h)(const double, const double),
		   const double x, const double y, double *out)
{
	impf_diff_1f2(f, x, y, out);
	impf_diff_1f2(g, x, y, out + 2);
	impf_diff_1f2(h, x, y, out + 4);
}
