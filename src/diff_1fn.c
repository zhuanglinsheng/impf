/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>
#include <string.h>

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

	memcpy(buffer, x, sizeof(double) * n);

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
