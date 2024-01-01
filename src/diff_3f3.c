/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>

void impf_diff_3f3(double (*f)(const double, const double, const double),
		   double (*g)(const double, const double, const double),
		   double (*h)(const double, const double, const double),
		   const double x, const double y, const double z, double *out)
{
	impf_diff_1f3(f, x, y, z, out);
	impf_diff_1f3(g, x, y, z, out + 3);
	impf_diff_1f3(h, x, y, z, out + 6);
}
