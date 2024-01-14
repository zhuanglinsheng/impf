/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

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
