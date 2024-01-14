/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#ifndef __IMPF_LINALG_H__
#define __IMPF_LINALG_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/*******************************************************************************
 * Linear algebra of "pivot-family"
 ******************************************************************************/

void impf_linalg_dscal(const int n, const double x, double *arr, const int inc);
void impf_linalg_daxpy(const int n, const double a, const double *x, const int incx, double *y, const int incy);
void impf_linalg_ldscal(const int n, const long double x, long double *arr, const int inc);
void impf_linalg_ldaxpy(const int n, const long double a, const long double *x, const int incx,
			long double *y, const int incy);

/*******************************************************************************
 * Linear algebra of "dgesl-family"
 ******************************************************************************/

int impf_linalg_dgels_2(const double A11, const double A12, const double A21, const double A22,
			const double b1, const double b2, double * x, int *code);
int impf_linalg_dgels_3(const double A11, const double A12, const double A13,
			const double A21, const double A22, const double A23,
			const double A31, const double A32, const double A33,
			const double b1, const double b2, const double b3,
			double *x, int *code);
/* Note: `buffer` is a double array, length = n * (n + nrhs) */
int impf_linalg_dgels_n(const int n, const int nrhs, const double *a, const int lda, const double *b,
			const int ldb, int *ipvt, double *x, double *buffer, int *code);
/* wrapper of `impf_linalg_dgels_n` by using impf_malloc */
int impf_linalg_dgels_L(const int n, const int nrhs, const double *a, const int lda, const double *b,
			const int ldb, int *ipvt, double *x, int *code);

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
