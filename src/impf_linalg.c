/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <assert.h>
#include <float.h>
#include <impf/impf_basic_blas.h>
#include <impf/impf_linalg.h>

impf_t_subrtstate impf_subrt_lsolve_2A2_b2(
    const double A11, const double A12, const double A21, const double A22,
    const double b1, const double b2, double * x)
{
    assert(NULL != x);

    double det = A11 * A22 - A12 * A21;
    *(x    ) = (A22 * b1 - A12 * b2) / det;
    *(x + 1) = (A11 * b2 - A21 * b1) / det;

    if(0 == det)
    {
        return IMPF_EXIT_FAILURE;
    }
    else if(impf_df_abs(det) < DBL_EPSILON)
    {
        return IMPF_EXIT_ILLNESS;
    }
    else
    {
        return IMPF_EXIT_SUCCESS;
    }
}

impf_t_subrtstate impf_subrt_lsolve_3A3_b3(
    const double A11, const double A12, const double A13,
    const double A21, const double A22, const double A23,
    const double A31, const double A32, const double A33,
    const double b1, const double b2, const double b3, double * x)
{
    assert(NULL != x);

    double detA, detA1, detA2, detA3;
    double tmp1, tmp2, tmp3;

    tmp1 = A22 * A33 - A23 * A32;
    tmp2 = A21 * A33 - A23 * A31;
    tmp3 = A21 * A32 - A22 * A31;
    detA = A11 * tmp1 - A12 * tmp2 + A13 * tmp3;

    detA1 = b1 * tmp1 - A12 * (b2 * A33 - b3 * A23) + A13 * (b2 * A32 - b3 * A22);
    detA2 = A11 * (b2 * A33 - b3 * A23) - b1 * tmp2 + A13 * (A21 * b3 - A31 * b2);
    detA3 = A11 * (A22 * b3 - A32 * b2) - A12 * (A21 * b3 - A31 * b2) + b1 * tmp3;

    *(x    ) = detA1 / detA;
    *(x + 1) = detA2 / detA;
    *(x + 2) = detA3 / detA;

    if(0 == detA)
    {
        return IMPF_EXIT_FAILURE;
    }
    else if(impf_df_abs(detA) < DBL_EPSILON)
    {
        return IMPF_EXIT_ILLNESS;
    }
    else
    {
        return IMPF_EXIT_SUCCESS;
    }
}
