/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/root.h>

int impf_root_1f1_secant(double (*f)(const double), double *x, double *y, int *code,
			 const double tol, const double tolf, const int niter)
{
	double fx, fy, z;
	int epoch;

	assert(f != NULL);
	assert(x != NULL);
	assert(y != NULL);
	assert(code != NULL);
	assert(tol  > 0);
	assert(tolf > 0);

	for (epoch = 0; epoch <= niter; epoch++) {
		fx = f(*x);
		fy = f(*y);

		if (tol > __impf_ABS__(*x - *y) && tolf > __impf_ABS__(fx)) {
			*code = impf_Success;
			return impf_EXIT_SUCCESS;
		}
		z = *y - ((*y) - (*x)) * fy / (fy - fx);
		*x = *y;
		*y = z;
	}
	*code = impf_ExceedIterLimit;
	return impf_EXIT_FAILURE;
}
