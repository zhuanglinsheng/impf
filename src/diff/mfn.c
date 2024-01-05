/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>
#include <string.h>

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

	memcpy(point, x, sizeof(double) * n);

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
