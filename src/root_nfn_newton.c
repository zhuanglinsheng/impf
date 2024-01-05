/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>
#include <impf/linalg.h>
#include <impf/root.h>
#include <string.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

static void diff_1f1(void (*f)(const double*, double*), double *x, double * Jac)
{
	double p1, p2, p3, p4;
	double f1, f2, f3, f4;
	double h = __impf_DIFF_H__;
	double two_h = 2. * h;

	p1 = *x - two_h;
	p2 = *x - h;
	p3 = *x + h;
	p4 = *x + two_h;
	f(&p1, &f1); f(&p2, &f2); f(&p3, &f3); f(&p4, &f4);
	Jac[0] = (f1 - f4 + 8. * (f3 - f2)) / (12. * h);
}

static void diff_2f2(void (*f)(const double*, double*), double *x, double * Jac)
{
	double p1[2], p2[2], p3[2], p4[2];
	double f1[2], f2[2], f3[2], f4[2];
	double h = __impf_DIFF_H__;
	double two_h = 2. * h;
	double h12 = 12. * h;

	p1[0] = x[0] - two_h; p1[1] = x[1];
	p2[0] = x[0] - h;     p2[1] = x[1];
	p3[0] = x[0] + h;     p3[1] = x[1];
	p4[0] = x[0] + two_h; p4[1] = x[1];
	f(p1, f1); f(p2, f2); f(p3, f3); f(p4, f4);
	Jac[0] = (f1[0] - f4[0] + 8. * (f3[0] - f2[0])) / h12;
	Jac[2] = (f1[1] - f4[1] + 8. * (f3[1] - f2[1])) / h12;
	p1[0] = x[0]; p1[1] = x[1] - two_h;
	p2[0] = x[0]; p2[1] = x[1] - h;
	p3[0] = x[0]; p3[1] = x[1] + h;
	p4[0] = x[0]; p4[1] = x[1] + two_h;
	f(p1, f1); f(p2, f2); f(p3, f3); f(p4, f4);
	Jac[1] = (f1[0] - f4[0] + 8. * (f3[0] - f2[0])) / h12;
	Jac[3] = (f1[1] - f4[1] + 8. * (f3[1] - f2[1])) / h12;
}

static void diff_3f3(void (*f)(const double*, double*), double *x, double * Jac)
{
	double p1[3], p2[3], p3[3], p4[3];
	double f1[3], f2[3], f3[3], f4[3];
	double h = __impf_DIFF_H__;
	double two_h = 2. * h;
	double h12 = 12. * h;

	p1[0] = x[0] - two_h; p1[1] = x[1]; p1[2] = x[2];
	p2[0] = x[0] - h;     p2[1] = x[1]; p2[2] = x[2];
	p3[0] = x[0] + h;     p3[1] = x[1]; p3[2] = x[2];
	p4[0] = x[0] + two_h; p4[1] = x[1]; p4[2] = x[2];
	f(p1, f1); f(p2, f2); f(p3, f3); f(p4, f4);
	Jac[0] = (f1[0] - f4[0] + 8. * (f3[0] - f2[0])) / h12;
	Jac[3] = (f1[1] - f4[1] + 8. * (f3[1] - f2[1])) / h12;
	Jac[6] = (f1[2] - f4[2] + 8. * (f3[2] - f2[2])) / h12;
	p1[0] = x[0]; p1[1] = x[1] - two_h;
	p2[0] = x[0]; p2[1] = x[1] - h;
	p3[0] = x[0]; p3[1] = x[1] + h;
	p4[0] = x[0]; p4[1] = x[1] + two_h;
	f(p1, f1); f(p2, f2); f(p3, f3); f(p4, f4);
	Jac[1] = (f1[0] - f4[0] + 8. * (f3[0] - f2[0])) / h12;
	Jac[4] = (f1[1] - f4[1] + 8. * (f3[1] - f2[1])) / h12;
	Jac[7] = (f1[2] - f4[2] + 8. * (f3[2] - f2[2])) / h12;
	p1[1] = x[1]; p1[2] = x[2] - two_h;
	p2[1] = x[1]; p2[2] = x[2] - h;
	p3[1] = x[1]; p3[2] = x[2] + h;
	p4[1] = x[1]; p4[2] = x[2] + two_h;
	f(p1, f1); f(p2, f2); f(p3, f3); f(p4, f4);
	Jac[2] = (f1[0] - f4[0] + 8. * (f3[0] - f2[0])) / h12;
	Jac[5] = (f1[1] - f4[1] + 8. * (f3[1] - f2[1])) / h12;
	Jac[8] = (f1[2] - f4[2] + 8. * (f3[2] - f2[2])) / h12;
}

int impf_root_nfn_newton(void (*f)(const double*, double*), double *x, const int n,
			 double *buffer, int *ipvt, int *code,
			 const double tol, const double tolf, const int niter)
{
	double *Jac, *x1_arr, *fx, *calc_buffer;
	int epoch;
	int i;

	assert(f != NULL);
	assert(x != NULL);
	assert(buffer != NULL);
	assert(code != NULL);
	assert(tol > 0);
	assert(tolf > 0);

	/*  buffer size = n(2n + 3). Deal with the cases n = 1, 2, 3 manually */
	Jac         = buffer;           /*       len = n * n                  */
	x1_arr      = Jac + n * n;      /*       len = n                      */
	fx          = x1_arr + n;       /*       len = n                      */
	calc_buffer = fx + n;           /*       len = n * max(    5,  n + 1) */
	f(x, fx);                       /* total len = n * max(n + 7, 2n + 3) */

	for (epoch = 0; epoch < niter; epoch++) {
		if (n == 1)
			diff_1f1(f, x, Jac);
		else if (n == 2)
			diff_2f2(f, x, Jac);
		else if (n == 3)
			diff_3f3(f, x, Jac);
		else
			impf_diff_mfn(f, x, n, Jac, n, calc_buffer);

		if (impf_linalg_dgels_n(n, 1, Jac, n, fx, 1, ipvt, x1_arr, calc_buffer, code) == EXIT_FAILURE)
			return EXIT_FAILURE; /* error code already updated*/
		for (i = 0; i < n; i++) {
			x1_arr[i] -= x[i];
			x1_arr[i] = -x1_arr[i];
		}
		f(x1_arr, fx);
		if (tolf >= maxabs_arrd(fx, n, 1)
		 && tol >= maxabs_arrd_gap(x, x1_arr, n, 1)) {
			memcpy(x, x1_arr, n * sizeof(double));
			*code = impf_Success;
			return EXIT_SUCCESS;
		}
		memcpy(x, x1_arr, n * sizeof(double));
	}
	*code = impf_ExceedIterLimit;
	return EXIT_FAILURE;
}
