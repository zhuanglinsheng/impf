/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>
#include <impf/root.h>

int impf_root_1f1_newton(double (*f)(const double), double *x, int *code,
			 const double tol, const double tolf, const int niter)
{
	double x1, grd, fx;
	int epoch;

	assert(f != NULL);
	assert(x != NULL);
	assert(code != NULL);
	assert(tol  > 0);
	assert(tolf > 0);

	fx = f(*x);

	for (epoch = 0; epoch <= niter; epoch++) {
		impf_diff_1f1(f, *x, &grd);

		if (__impf_IDF_LINALG_DET_ZERO__ >= __impf_ABS__(grd)) {
			*code = impf_Singularity;
			return impf_EXIT_FAILURE;
		}
		x1 = *x - fx / grd;
		fx = f(x1);

		if (tol >= __impf_ABS__(*x - x1) && tolf >= __impf_ABS__(fx)) {
			*x = x1;
			*code = impf_Success;
			return impf_EXIT_SUCCESS;
		}
		*x = x1;
	}
	*code = impf_ExceedIterLimit;
	return impf_EXIT_FAILURE;
}
