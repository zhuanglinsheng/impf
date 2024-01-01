/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>

int impf_linalg_dgels_2(const double A11, const double A12, const double A21, const double A22,
			const double b1, const double b2, double *x, enum impf_ErrorCode *code)
{
	double det;

	assert(x != NULL);

	det = A11 * A22 - A12 * A21;
	*(x) = (A22 * b1 - A12 * b2) / det;
	*(x + 1) = (A11 * b2 - A21 * b1) / det;

	if (__impf_ABS__(det) < __impf_IDF_LINALG_DET_ZERO__) {
		*code = impf_Singularity;
		return EXIT_FAILURE;
	} else {
		*code = impf_Success;
		return EXIT_SUCCESS;
	}
}
