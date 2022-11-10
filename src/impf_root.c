/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <float.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <impf/impf_basic_rfs.h>
#include <impf/impf_basic.h>
#include <impf/impf_ndiff.h>
#include <impf/impf_linalg.h>
#include <impf/impf_root.h>

impf_t_subrtstate
impf_subrt_root_newton_1f1(
    impf_t_f1 f, double * x, const double tol, const double tolf, const unsigned int niter)
{
    double x0, x1, grd, fx;
    unsigned int epoch;

    x0 = *x;
    fx = f(x0);
    epoch = 0;

LOOP:
    impf_subrt_diff(f, x0, &grd);

    if(DBL_EPSILON >= fabs(grd))
    {
        return IMPF_EXIT_FAILURE;
    }
    x1 = x0 - fx / grd;
    fx = f(x1);
    epoch ++;

    if(epoch >= niter)
    {
        *x = x1;
        return IMPF_EXIT_ILLNESS;
    }
    if(tol >= impf_df_abs_sub(x0, x1) && tolf >= fabs(fx))
    {
        *x = x1;
        return IMPF_EXIT_SUCCESS;
    }
    x0 = x1;
    goto LOOP;
}

impf_t_subrtstate
impf_subrt_root_newton_2f2(
    impf_t_f2 f1, impf_t_f2 f2,
    double * x, double * y,
    const double tol, const double tolf, const unsigned int niter)
{
    double jac_data[2 * 2];
    double x0_array[2];
    double x1_array[2];
    double f12_xy[2];
    impf_t_matrix jac = {jac_data, 2, 2, IMPF_MAT_ROW_MAJOR};

    unsigned int epoch;

    x0_array[0] = *x;
    x0_array[1] = *y;
    f12_xy[0] = f1(x0_array[0], x0_array[1]);
    f12_xy[1] = f2(x0_array[0], x0_array[1]);
    epoch = 0;

LOOP:
    /* Calculate J */
    impf_subrt_jacob_2f2(f1, f2, x0_array[0], x0_array[1], &jac);

    /* Solve J^{-1} f12_xy^T */
    if(IMPF_EXIT_SUCCESS != impf_subrt_lsolve_2A2_b2(
        jac_data[0], jac_data[1], \
        jac_data[2], jac_data[3], \
        f12_xy[0], f12_xy[1], x1_array))
    {
        return IMPF_EXIT_FAILURE;
    }

    /* Newton update of solution */
    x1_array[0] = x0_array[0] - x1_array[0];
    x1_array[1] = x0_array[1] - x1_array[1];

    /* Updated values and epochs */
    f12_xy[0] = f1(x1_array[0], x1_array[1]);
    f12_xy[1] = f2(x1_array[0], x1_array[1]);
    epoch ++;

    if(epoch >= niter)
    {
        *x = x1_array[0];
        *y = x1_array[1];
        return IMPF_EXIT_ILLNESS;
    }
    if(tol >= impf_df_distance_norminf_2d(x0_array, x1_array)
    && tolf >= impf_df_norminf_2d(f12_xy))
    {
        *x = x1_array[0];
        *y = x1_array[1];
        return IMPF_EXIT_SUCCESS;
    }
    x0_array[0] = x1_array[0];
    x0_array[1] = x1_array[1];
    goto LOOP;
}

impf_t_subrtstate
impf_subrt_root_newton_3f3(
    impf_t_f3 f1, impf_t_f3 f2, impf_t_f3 f3,
    double * x, double * y, double * z,
    const double tol, const double tolf, const unsigned int niter)
{
    double jac_data[3 * 3];
    double x0_array[3];
    double x1_array[3];
    double f123_xyz[3];
    impf_t_matrix jac = {jac_data, 3, 3, IMPF_MAT_ROW_MAJOR};

    unsigned int epoch;

    x0_array[0] = *x;
    x0_array[1] = *y;
    x0_array[2] = *z;
    f123_xyz[0] = f1(x0_array[0], x0_array[1], x0_array[2]);
    f123_xyz[1] = f2(x0_array[0], x0_array[1], x0_array[2]);
    f123_xyz[2] = f3(x0_array[0], x0_array[1], x0_array[2]);
    epoch = 0;

LOOP:
    /* Calculate J */
    impf_subrt_jacob_3f3(f1, f2, f3, x0_array[0], x0_array[1], x0_array[2], &jac);

    /* Solve J^{-1} f123_xyz^T */
    if(IMPF_EXIT_SUCCESS != impf_subrt_lsolve_3A3_b3(
        jac_data[0], jac_data[1], jac_data[2], \
        jac_data[3], jac_data[4], jac_data[5], \
        jac_data[6], jac_data[7], jac_data[8], \
        f123_xyz[0], f123_xyz[1], f123_xyz[2], x1_array))
    {
        return IMPF_EXIT_FAILURE;
    }

    /* Newton update of solution */
    x1_array[0] = x0_array[0] - x1_array[0];
    x1_array[1] = x0_array[1] - x1_array[1];
    x1_array[2] = x0_array[2] - x1_array[2];

    /* Updated values and epochs */
    f123_xyz[0] = f1(x1_array[0], x1_array[1], x1_array[2]);
    f123_xyz[1] = f2(x1_array[0], x1_array[1], x1_array[2]);
    f123_xyz[2] = f3(x1_array[0], x1_array[1], x1_array[2]);
    epoch ++;

    if(epoch >= niter)
    {
        *x = x1_array[0];
        *y = x1_array[1];
        *z = x1_array[2];
        return IMPF_EXIT_ILLNESS;
    }
    if(tol >= impf_df_distance_norminf_3d(x0_array, x1_array)
    && tolf >= impf_df_norminf_3d(f123_xyz))
    {
        *x = x1_array[0];
        *y = x1_array[1];
        *z = x1_array[2];
        return IMPF_EXIT_SUCCESS;
    }
    x0_array[0] = x1_array[0];
    x0_array[1] = x1_array[1];
    x0_array[2] = x1_array[2];
    goto LOOP;
}
