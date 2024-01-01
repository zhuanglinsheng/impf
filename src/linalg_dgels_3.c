/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>

int impf_linalg_dgels_3(const double A11, const double A12, const double A13,
			const double A21, const double A22, const double A23,
			const double A31, const double A32, const double A33,
			const double b1, const double b2, const double b3, double *x, enum impf_ErrorCode *code)
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
		return EXIT_FAILURE;
	}
	detA1 = b1 * tmp1 - A12 * (b2 * A33 - b3 * A23) + A13 * (b2 * A32 - b3 * A22);
	detA2 = A11 * (b2 * A33 - b3 * A23) - b1 * tmp2 + A13 * (A21 * b3 - A31 * b2);
	detA3 = A11 * (A22 * b3 - A32 * b2) - A12 * (A21 * b3 - A31 * b2) + b1 * tmp3;

	*(x)     = detA1 / detA;
	*(x + 1) = detA2 / detA;
	*(x + 2) = detA3 / detA;
	*code = impf_Success;
	return EXIT_SUCCESS;
}
