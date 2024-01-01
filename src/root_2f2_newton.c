/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>
#include <impf/root.h>
#include <impf/linalg.h>

int impf_root_2f2_newton(double (*f)(const double, const double), double (*g)(const double, const double),
			 double *x, double *y, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter)
{
	double jac[4], x1_arr[2];
	double f0, g0;
	int epoch;

	assert(f != NULL);
	assert(g != NULL);
	assert(x != NULL);
	assert(y != NULL);
	assert(code != NULL);
	assert(tol  > 0);
	assert(tolf > 0);

	f0 = f(*x, *y);
	g0 = g(*x, *y);

	for (epoch = 0; epoch <= niter; epoch++) {
		impf_diff_2f2(f, g, *x, *y, jac);

		if (impf_linalg_dgels_2(jac[0], jac[1], jac[2], jac[3], f0, g0, x1_arr, code) == EXIT_FAILURE)
			return EXIT_FAILURE; /* code already updated*/
		x1_arr[0] = *x - x1_arr[0];
		x1_arr[1] = *y - x1_arr[1];
		f0 = f(x1_arr[0], x1_arr[1]);
		g0 = g(x1_arr[0], x1_arr[1]);

		if (tol  >= __impf_ABS__(*x - x1_arr[0])
		 && tol  >= __impf_ABS__(*y - x1_arr[1])
		 && tolf >= __impf_ABS__(f0)
		 && tolf >= __impf_ABS__(g0)) {
			*x = x1_arr[0];
			*y = x1_arr[1];
			*code = impf_Success;
			return EXIT_SUCCESS;
		}
		*x = x1_arr[0];
		*y = x1_arr[1];
	}
	*code = impf_ExceedIterLimit;
	return EXIT_FAILURE;
}
