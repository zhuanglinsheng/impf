/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>

void impf_diff_1f3(double (*f)(const double, const double, const double),
		   const double x, const double y, const double z, double *out)
{
	double h, two_h, twelf_h;
	double x1, x2, x3, x4;
	double y1, y2, y3, y4;
	double z1, z2, z3, z4;

	assert(f != NULL);
	assert(out != NULL);

	h = __impf_DIFF_H__;
	two_h = 2. * h;
	twelf_h = 12. * h;

	x1 = x - two_h;
	y1 = y - two_h;
	z1 = z - two_h;
	x2 = x - h;
	y2 = y - h;
	z2 = z - h;
	x3 = x + h;
	y3 = y + h;
	z3 = z + h;
	x4 = x + two_h;
	y4 = y + two_h;
	z4 = z + two_h;

	*out       = (f(x1, y, z) - f(x4, y, z) + 8. * (f(x3, y, z) - f(x2, y, z))) / twelf_h;
	*(out + 1) = (f(x, y1, z) - f(x, y4, z) + 8. * (f(x, y3, z) - f(x, y2, z))) / twelf_h;
	*(out + 2) = (f(x, y, z1) - f(x, y, z4) + 8. * (f(x, y, z3) - f(x, y, z2))) / twelf_h;
}
