/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#ifndef __IMPF_DIFF_H__
#define __IMPF_DIFF_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/*******************************************************************************
 * Numerical differentiation of "diff-first-order-family"
 *
 * Assuption:
 *	1. for any given x, f(x) is accurate
 *	2. points selection (p1, p2, p3, p4) are accurate (x is NOT too large)
 *
 * Theoretical error (five-point formula) = h^4 * f^(5)(x) / 30
 * Float rounding error = 3 eps / (4 * h)
 * Optimal h is approximately eps^{1/5}
 ******************************************************************************/

void impf_diff_1f1(double (*f)(const double), const double x, double *out);
void impf_diff_2f1(double (*f)(const double),
		   double (*g)(const double), const double x, double *out);
void impf_diff_3f1(double (*f)(const double), double (*g)(const double),
		   double (*h)(const double), const double x, double *out);

void impf_diff_1f2(double (*f)(const double, const double),
		   const double x, const double y, double *out);
void impf_diff_2f2(double (*f)(const double, const double),
		   double (*g)(const double, const double),
		   const double x, const double y, double *out);
void impf_diff_3f2(double (*f)(const double, const double),
		   double (*g)(const double, const double),
		   double (*h)(const double, const double),
		   const double x, const double y, double *out);

void impf_diff_1f3(double (*f)(const double, const double, const double),
		   const double x, const double y, const double z, double *out);
void impf_diff_2f3(double (*f)(const double, const double, const double),
		   double (*g)(const double, const double, const double),
		   const double x, const double y, const double z, double *out);
void impf_diff_3f3(double (*f)(const double, const double, const double),
		   double (*g)(const double, const double, const double),
		   double (*h)(const double, const double, const double),
		   const double x, const double y, const double z, double *out);

/* Note: `buffer` is a double array, length = n */
void impf_diff_1fn(double (*f)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer);
/* Note: `buffer` is a double array, length = n */
void impf_diff_2fn(double (*f)(const double*, const int),
		   double (*g)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer);
/* Note: `buffer` is a double array, length = n */
void impf_diff_3fn(double (*f)(const double*, const int),
		   double (*g)(const double*, const int),
		   double (*h)(const double*, const int),
		   const double *x, const int n, double *out, double *buffer);
/* Note: `buffer` is a double array, length = n + 4m */
void impf_diff_mfn(void (*f)(const double*, double*), const double *x, const int n,
		   double *out, const int m, double *buffer);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
