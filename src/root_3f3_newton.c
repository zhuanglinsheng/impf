/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/diff.h>
#include <impf/root.h>
#include <impf/linalg.h>

int impf_root_3f3_newton(double (*f)(const double, const double, const double),
			 double (*g)(const double, const double, const double),
			 double (*h)(const double, const double, const double),
			 double *x, double *y, double *z, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter)
{
	double jac[9], x1_arr[3];
	double f0, g0, h0;
	int epoch;

	assert(f != NULL);
	assert(g != NULL);
	assert(h != NULL);
	assert(x != NULL);
	assert(y != NULL);
	assert(z != NULL);
	assert(code != NULL);
	assert(tol > 0);
	assert(tolf > 0);

	f0 = f(*x, *y, *z);
	g0 = g(*x, *y, *z);
	h0 = h(*x, *y, *z);

	for (epoch = 0; epoch < niter; epoch++) {
		impf_diff_3f3(f, g, h, *x, *y, *z, jac);

		if (impf_linalg_dgels_3(jac[0], jac[1], jac[2],
					jac[3], jac[4], jac[5],
					jac[6], jac[7], jac[8],
					f0, g0, h0, x1_arr, code) == EXIT_FAILURE)
			return EXIT_FAILURE; /* code already updated*/
		x1_arr[0] = *x - x1_arr[0];
		x1_arr[1] = *y - x1_arr[1];
		x1_arr[2] = *z - x1_arr[2];
		f0 = f(x1_arr[0], x1_arr[1], x1_arr[2]);
		g0 = g(x1_arr[0], x1_arr[1], x1_arr[2]);
		h0 = h(x1_arr[0], x1_arr[1], x1_arr[2]);

		if (tol  >= __impf_ABS__(*x - x1_arr[0])
		 && tol  >= __impf_ABS__(*y - x1_arr[1])
		 && tol  >= __impf_ABS__(*z - x1_arr[2])
		 && tolf >= __impf_ABS__(f0)
		 && tolf >= __impf_ABS__(g0)
		 && tolf >= __impf_ABS__(h0)) {
			*x = x1_arr[0];
			*y = x1_arr[1];
			*z = x1_arr[2];
			*code = impf_Success;
			return EXIT_SUCCESS;
		}
		*x = x1_arr[0];
		*y = x1_arr[1];
		*z = x1_arr[2];
	}
	*code = impf_ExceedIterLimit;
	return EXIT_FAILURE;
}
