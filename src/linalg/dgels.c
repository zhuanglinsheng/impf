/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

int impf_linalg_dgels_2(const double A11, const double A12, const double A21, const double A22,
			const double b1, const double b2, double *x, int *code)
{
	double det;

	assert(x != NULL);

	det = A11 * A22 - A12 * A21;
	*(x) = (A22 * b1 - A12 * b2) / det;
	*(x + 1) = (A11 * b2 - A21 * b1) / det;

	if (__impf_ABS__(det) < __impf_IDF_LINALG_DET_ZERO__) {
		*code = impf_Singularity;
		return impf_EXIT_FAILURE;
	} else {
		*code = impf_Success;
		return impf_EXIT_SUCCESS;
	}
}

int impf_linalg_dgels_3(const double A11, const double A12, const double A13,
			const double A21, const double A22, const double A23,
			const double A31, const double A32, const double A33,
			const double b1, const double b2, const double b3, double *x, int *code)
{
	double detA, detA1, detA2, detA3;
	double tmp1, tmp2, tmp3;

	assert(x != NULL);
	assert(code != NULL);

	tmp1 = A22 * A33 - A23 * A32;
	tmp2 = A21 * A33 - A23 * A31;
	tmp3 = A21 * A32 - A22 * A31;
	detA = A11 * tmp1 - A12 * tmp2 + A13 * tmp3;

	if (__impf_ABS__(detA) < __impf_IDF_LINALG_DET_ZERO__) {
		*code = impf_Singularity;
		return impf_EXIT_FAILURE;
	}
	detA1 = b1 * tmp1 - A12 * (b2 * A33 - b3 * A23) + A13 * (b2 * A32 - b3 * A22);
	detA2 = A11 * (b2 * A33 - b3 * A23) - b1 * tmp2 + A13 * (A21 * b3 - A31 * b2);
	detA3 = A11 * (A22 * b3 - A32 * b2) - A12 * (A21 * b3 - A31 * b2) + b1 * tmp3;

	*(x)     = detA1 / detA;
	*(x + 1) = detA2 / detA;
	*(x + 2) = detA3 / detA;
	*code = impf_Success;
	return impf_EXIT_SUCCESS;
}

int impf_linalg_dgels_L(const int n, const int nrhs, const double *a, const int lda, const double *b, const int ldb,
			int *ipvt, double *x, int *code)
{
	double *buffer = impf_malloc(n * (n + nrhs) * sizeof(double));

	if (buffer == NULL) {
		*code = impf_MemoryAllocError;
		return impf_EXIT_FAILURE;
	}
	return impf_linalg_dgels_n(n, nrhs, a, lda, b, ldb, ipvt, x, buffer, code);
}

static int select_major(const double *row, const int n, const int *ipvt, const int i)
{
	int j, idx = 0;
	double ele, maxv = __impf_NINF__;

	for (j = 0; j < n; j++) {
		if (is_in_arri(j, ipvt, i))
			continue;
		ele = __impf_ABS__(row[j]);
		if (ele > maxv) {
			maxv = ele;
			idx = j;
		}
	}
	return idx;
}

int impf_linalg_dgels_n(const int n, const int nrhs, const double *a, const int lda, const double *b, const int ldb,
			int *ipvt, double *x, double *buffer, int *code)
{
	int i, j;
	int ncol = n + nrhs;

	assert(a != NULL);
	assert(ipvt != NULL);
	assert(b != NULL);
	assert(x != NULL);
	assert(ipvt != NULL);
	assert(buffer != NULL);
	assert(code != NULL);

	for (i = 0; i < n; i++) {
		int rowi = i * ncol;

		impf_memcpy(buffer + rowi, a + i * lda, n * sizeof(double));
		impf_memcpy(buffer + rowi + n, b + i * ldb, nrhs * sizeof(double));
	}
	impf_memset(ipvt, 0, n * sizeof(int));

	for (i = 0; i < n; i++) {
		int rowi = i * ncol;
		double *prowi = buffer + rowi;
		int majorcol = select_major(prowi, n, ipvt, i);
		double majorv_i = buffer[majorcol + rowi];

		ipvt[i] = majorcol;

		if (__impf_ABS__(majorv_i) < __impf_IDF_LINALG_MAJ_ZERO__) {
			*code = impf_Singularity;
			return impf_EXIT_FAILURE;
		}
		impf_linalg_dscal(ncol, 1. / majorv_i, prowi, 1);

		for (j = 0; j < n; j++) {
			double *prowj = buffer + j * ncol;
			double majorv_j = buffer[majorcol + j * ncol];

			if (j == i)
				continue;
			impf_linalg_daxpy(ncol, -majorv_j, prowi, 1, prowj, 1);
		}
	}
	for (i = 0; i < n; i++) {
		int rowi = i * ncol;
		int rowx = ipvt[i] * nrhs;

		impf_memcpy(x + rowx, buffer + n + rowi, nrhs * sizeof(double));
	}
	*code = impf_Success;
	return impf_EXIT_SUCCESS;
}
