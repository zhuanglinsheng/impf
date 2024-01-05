/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>

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
