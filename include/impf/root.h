/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#ifndef __IMPF_ROOT_H__
#define __IMPF_ROOT_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/*******************************************************************************
 * Root finding of "root-1f1-family"
 ******************************************************************************/

/* Root finding: 1-variant real function (bisection method)
 *
 * Note:
 * 	1. `x` and `y` are initial boundaries of solution
 *	2. the initial boundaries `x` and `y` shall satisfy f(x) * f(y) < 0
 */
int impf_root_1f1_bisection(double (*f)(const double), double *x, double *y, enum impf_ErrorCode *code,
			    const double tol, const double tolf, const int niter);

/* Root finding: 1-variant real function (secant method)
 *
 * Note: `x` and `y` are initial guesses of solution
 */
int impf_root_1f1_secant(double (*f)(const double), double *x, double *y, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter);

/* Root finding: 1-variant real function (newton's method)
 *
 * Note: `x` is the initial guess of solution
 */
int impf_root_1f1_newton(double (*f)(const double), double *x, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter);

/*******************************************************************************
 * Root finding of "root-nfn-family"
 ******************************************************************************/

/* Root finding: two 2-variant real functions (newton's method)
 *
 * Note: `x` and `y` are the initial guesses of solution
 */
int impf_root_2f2_newton(double (*f)(const double, const double),
			 double (*g)(const double, const double),
			 double *x, double *y, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter);

/* Root finding: three 3-variant real functions (newton's method)
 *
 * Note: `x`, `y` and `z` are the initial guesses of solution
 */
int impf_root_3f3_newton(double (*f)(const double, const double, const double),
			 double (*g)(const double, const double, const double),
			 double (*h)(const double, const double, const double),
			 double *x, double *y, double *z, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter);

/* Root finding: n-variant real function that returns n values (newton's method)
 *
 * Note: (1) `x` is a double array providing the initial guess of solution, length = n
 *       (2) `buffer` is a double array, length = n(2n + 3)
 *       (3) `pivt` is an int array, length = n
 */
int impf_root_nfn_newton(void (*f)(const double*, double*), double *x, const int n,
			 double *buffer, int *pivt, enum impf_ErrorCode *code,
			 const double tol, const double tolf, const int niter);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
