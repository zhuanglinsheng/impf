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
