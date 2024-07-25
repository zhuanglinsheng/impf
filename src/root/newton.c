/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/diff.h>
#include <impf/linalg.h>
#include <impf/root.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

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

int impf_root_2f2_newton(double (*f)(const double, const double), double (*g)(const double, const double),
			 double *x, double *y, int *code,
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

		if (impf_linalg_dgels_2(jac[0], jac[1], jac[2], jac[3], f0, g0, x1_arr, code) == impf_EXIT_FAILURE)
			return impf_EXIT_FAILURE; /* code already updated*/
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
			return impf_EXIT_SUCCESS;
		}
		*x = x1_arr[0];
		*y = x1_arr[1];
	}
	*code = impf_ExceedIterLimit;
	return impf_EXIT_FAILURE;
}

int impf_root_3f3_newton(double (*f)(const double, const double, const double),
			 double (*g)(const double, const double, const double),
			 double (*h)(const double, const double, const double),
			 double *x, double *y, double *z, int *code,
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
					f0, g0, h0, x1_arr, code) == impf_EXIT_FAILURE)
			return impf_EXIT_FAILURE; /* code already updated*/
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
			return impf_EXIT_SUCCESS;
		}
		*x = x1_arr[0];
		*y = x1_arr[1];
		*z = x1_arr[2];
	}
	*code = impf_ExceedIterLimit;
	return impf_EXIT_FAILURE;
}

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

		if (impf_linalg_dgels_n(n, 1, Jac, n, fx, 1, ipvt, x1_arr, calc_buffer, code) == impf_EXIT_FAILURE)
			return impf_EXIT_FAILURE; /* error code already updated*/
		for (i = 0; i < n; i++) {
			x1_arr[i] -= x[i];
			x1_arr[i] = -x1_arr[i];
		}
		f(x1_arr, fx);
		if (tolf >= maxabs_arrd(fx, n, 1)
		 && tol >= maxabs_arrd_gap(x, x1_arr, n, 1)) {
			impf_memcpy(x, x1_arr, n * sizeof(double));
			*code = impf_Success;
			return impf_EXIT_SUCCESS;
		}
		impf_memcpy(x, x1_arr, n * sizeof(double));
	}
	*code = impf_ExceedIterLimit;
	return impf_EXIT_FAILURE;
}
