/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/linalg.h>

int impf_linalg_dgels_L(const int n, const int nrhs, const double *a, const int lda, const double *b, const int ldb,
			int *ipvt, double *x, int *code)
{
	double *buffer = malloc(n * (n + nrhs) * sizeof(double));

	if (buffer == NULL) {
		*code = impf_MemoryAllocError;
		return EXIT_FAILURE;
	}
	return impf_linalg_dgels_n(n, nrhs, a, lda, b, ldb, ipvt, x, buffer, code);
}
