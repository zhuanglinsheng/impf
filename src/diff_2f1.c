/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>

void impf_diff_2f1(double (*f)(const double), double (*g)(const double), const double x, double *out)
{
	impf_diff_1f1(f, x, out);
	impf_diff_1f1(g, x, out + 1);
}
