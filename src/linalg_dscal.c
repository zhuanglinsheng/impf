/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/linalg.h>

void impf_linalg_dscal(const int n, const double x, double *arr, const int inc)
{
	int i;

	assert(arr != NULL);

	for (i = 0; i < n; i += inc)
		arr[i] *= x;
}
