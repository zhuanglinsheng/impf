/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <float.h>
#include <string.h>
#include <assert.h>
#include <impf/impf_basic.h>
#include <impf/impf_ndiff.h>

/**
 * \brief The maximum DBL_EPSILON by ISO C90 standard.
 *
 * \details The C90 standard gives a maximum DBL_EPSILON in <float.h>, which is
 * not necessarily the smallest value that can be expressed by the machine.
 */
#define __impf_C_MAX_DBL_EPS__  1e-9

#define __impf_SQRT_DBL_EPS_V1__ ((__impf_C_MAX_DBL_EPS__ + DBL_EPSILON / __impf_C_MAX_DBL_EPS__)/2.)
#define __impf_SQRT_DBL_EPS_V2__ ((__impf_SQRT_DBL_EPS_V1__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V1__)/2.)
#define __impf_SQRT_DBL_EPS_V3__ ((__impf_SQRT_DBL_EPS_V2__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V2__)/2.)
#define __impf_SQRT_DBL_EPS_V4__ ((__impf_SQRT_DBL_EPS_V3__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V3__)/2.)
#define __impf_SQRT_DBL_EPS_V5__ ((__impf_SQRT_DBL_EPS_V4__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V4__)/2.)
#define __impf_SQRT_DBL_EPS_V6__ ((__impf_SQRT_DBL_EPS_V5__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V5__)/2.)
#define __impf_SQRT_DBL_EPS_V7__ ((__impf_SQRT_DBL_EPS_V6__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V6__)/2.)
#define __impf_SQRT_DBL_EPS_V8__ ((__impf_SQRT_DBL_EPS_V7__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V7__)/2.)
#define __impf_SQRT_DBL_EPS_V9__ ((__impf_SQRT_DBL_EPS_V8__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V8__)/2.)

/**
 * \brief The squared root of `DBL_EPSILON` defined in <float.h>.
 *
 * \details The constant is defined as a global variable, since it is machine
 * dependent. We use an ISO C90 standard defined maximum DBL_EPSILON as the
 * initial guess, and estimate the squared root using the Newton's iterations
 * ten times by
 * \f[
 * x_{t+1} = \frac{x_t+\epsilon/x_t}{2}
 * \f]
 * The macro `DBL_EPSILON` is defined in <float.h>, which is machine dependent.
 * The entity is defined in the file "impf_ndiff.c".
 */
const double impf_SQRT_DBL_EPSILON = \
    ((__impf_SQRT_DBL_EPS_V9__ + DBL_EPSILON / __impf_SQRT_DBL_EPS_V9__)/2.);

/* Cancels the definitions of middle steps */
#undef __impf_SQRT_DBL_EPS_V9__
#undef __impf_SQRT_DBL_EPS_V8__
#undef __impf_SQRT_DBL_EPS_V7__
#undef __impf_SQRT_DBL_EPS_V6__
#undef __impf_SQRT_DBL_EPS_V5__
#undef __impf_SQRT_DBL_EPS_V4__
#undef __impf_SQRT_DBL_EPS_V3__
#undef __impf_SQRT_DBL_EPS_V2__
#undef __impf_SQRT_DBL_EPS_V1__
#undef __impf_C_MAX_DBL_EPS__


void impf_subrt_diff(
    const impf_t_f1 f,
    const double x, 
    double * re)
{
    assert(NULL != f);
    assert(NULL != re);

    double h, two_h;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    two_h = 2 * h;
    *re = (f(x - two_h) - f(x + two_h) + 8 * (f(x + h) - f(x - h))) / (12 * h);
}


void impf_subrt_gradient_f2(
    const impf_t_f2 f,
    const double x, const double y,
    double * re)
{
    assert(NULL != f);
    assert(NULL != re);

    double h, two_h, twelf_h;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    two_h = 2 * h;
    twelf_h = 12 * h;
    *re       = (f(x - two_h, y) - f(x + two_h, y) + 8 * (f(x + h, y) - f(x - h, y))) / twelf_h;
    *(re + 1) = (f(x, y - two_h) - f(x, y + two_h) + 8 * (f(x, y + h) - f(x, y - h))) / twelf_h;
}


void impf_subrt_gradient_f3(
    const impf_t_f3 f,
    const double x, const double y, const double z,
    double * re)
{
    assert(NULL != f);
    assert(NULL != re);

    double h, two_h, twelf_h;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    two_h = 2 * h;
    twelf_h = 12 * h;
    *re       = (f(x - two_h, y, z) - f(x + two_h, y, z) + 8 * (f(x + h, y, z) - f(x - h, y, z))) / twelf_h;
    *(re + 1) = (f(x, y - two_h, z) - f(x, y + two_h, z) + 8 * (f(x, y + h, z) - f(x, y - h, z))) / twelf_h;
    *(re + 2) = (f(x, y, z - two_h) - f(x, y, z + two_h) + 8 * (f(x, y, z + h) - f(x, y, z - h))) / twelf_h;
}


void impf_subrt_gradient_fn(
    const impf_t_fn f,
    const double * x, const unsigned int n,
    double * buf, double * grd)
{
    assert(NULL != f);
    assert(NULL != x);
    assert(NULL != buf);
    assert(NULL != grd);

    double h, two_h, twelf_h, tmp1, tmp2, tmp3, tmp4;
    unsigned int i;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    two_h = 2 * h;
    twelf_h = 12 * h;
    memcpy(buf, x, sizeof(double) * n);  /* Copy x into buf */

    for(i = 0; i<n; i++)
    {
        *(buf + i) -= two_h;
        tmp1 = f(buf, n);     /* f(x-2h) */
        *(buf + i) += h;  
        tmp2 = f(buf, n);     /* f(x-h)  */
        *(buf + i) += two_h;
        tmp3 = f(buf, n);     /* f(x+h)  */
        *(buf + i) += h;
        tmp4 = f(buf, n);     /* f(x+2h) */
        *(buf + i) -= two_h;
        *(grd + i) = (tmp1 - tmp4 + 8 * (tmp3 - tmp2)) / twelf_h;
    }
}


void impf_subrt_gradient(
    const impf_t_fn f, const impf_t_vector * vec,
    double * buf, double * grd)
{
    assert(NULL != vec);
    impf_subrt_gradient_fn(f, vec->data, vec->dim, buf, grd);
}



void impf_subrt_jacob_2f1(
    const impf_t_f1 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f1 f2,      /* [in] the 2nd differentiated function. */ 
    const double x,          /* [in] the evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$2\times1\f$.  */
{
    assert(NULL != jac);

    double d1, d2;

    impf_subrt_diff(f1, x, &d1);
    impf_subrt_diff(f2, x, &d2);
    *(jac->data    ) = d1;
    *(jac->data + 1) = d2;
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 2;
    jac->ncol = 1;
}


void impf_subrt_jacob_2f2(
    const impf_t_f2 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f2 f2,      /* [in] the 2nd differentiated function. */ 
    const double x,          /* [in] the 1st evaluation point. */
    const double y,          /* [in] the 2nd evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$2\times2\f$. */
{
    assert(NULL != jac);

    double df1[2], df2[2];

    impf_subrt_gradient_f2(f1, x, y, df1);
    impf_subrt_gradient_f2(f2, x, y, df2);
    memcpy(jac->data    , df1, sizeof(double) * 2);
    memcpy(jac->data + 2, df2, sizeof(double) * 2);
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 2;
    jac->ncol = 2;
}


void impf_subrt_jacob_2f3(
    const impf_t_f3 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f3 f2,      /* [in] the 2nd differentiated function. */ 
    const double x,          /* [in] the 1st evaluation point. */
    const double y,          /* [in] the 2nd evaluation point. */
    const double z,          /* [in] the 3rd evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$2\times3\f$. */
{
    assert(NULL != jac);

    double df1[3], df2[3];

    impf_subrt_gradient_f3(f1, x, y, z, df1);
    impf_subrt_gradient_f3(f2, x, y, z, df2);
    memcpy(jac->data    , df1, sizeof(double) * 3);
    memcpy(jac->data + 3, df2, sizeof(double) * 3);
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 2;
    jac->ncol = 3;
}


void impf_subrt_jacob_3f1(
    const impf_t_f1 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f1 f2,      /* [in] the 2nd differentiated function. */ 
    const impf_t_f1 f3,      /* [in] the 3rd differentiated function. */ 
    const double x,          /* [in] the evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$3\times1\f$.  */
{
    assert(NULL != jac);

    double df1, df2, df3;

    impf_subrt_diff(f1, x, &df1);
    impf_subrt_diff(f2, x, &df2);
    impf_subrt_diff(f3, x, &df3);
    *(jac->data    ) = df1;
    *(jac->data + 1) = df2;
    *(jac->data + 2) = df3;
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 3;
    jac->ncol = 1;
}


void impf_subrt_jacob_3f2(
    const impf_t_f2 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f2 f2,      /* [in] the 2nd differentiated function. */ 
    const impf_t_f2 f3,      /* [in] the 3rd differentiated function. */ 
    const double x,          /* [in] the 1st evaluation point. */
    const double y,          /* [in] the 2nd evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$3\times2\f$. */
{
    assert(NULL != jac);

    double df1[2], df2[2], df3[2];

    impf_subrt_gradient_f2(f1, x, y, df1);
    impf_subrt_gradient_f2(f2, x, y, df2);
    impf_subrt_gradient_f2(f3, x, y, df3);
    memcpy(jac->data    , df1, sizeof(double) * 2);
    memcpy(jac->data + 2, df2, sizeof(double) * 2);
    memcpy(jac->data + 4, df3, sizeof(double) * 2);
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 3;
    jac->ncol = 2;
}


void impf_subrt_jacob_3f3(
    const impf_t_f3 f1,      /* [in] the 1st differentiated function. */ 
    const impf_t_f3 f2,      /* [in] the 2nd differentiated function. */ 
    const impf_t_f3 f3,      /* [in] the 3rd differentiated function. */ 
    const double x,          /* [in] the 1st evaluation point. */
    const double y,          /* [in] the 2nd evaluation point. */
    const double z,          /* [in] the 3rd evaluation point. */
    impf_t_matrix * jac)     /* [out] the Jacobian matrix of dimension \f$3\times3\f$. */
{
    assert(NULL != jac);

    double df1[3], df2[3], df3[3];

    impf_subrt_gradient_f3(f1, x, y, z, df1);
    impf_subrt_gradient_f3(f2, x, y, z, df2);
    impf_subrt_gradient_f3(f3, x, y, z, df3);
    memcpy(jac->data    , df1, sizeof(double) * 3);
    memcpy(jac->data + 3, df2, sizeof(double) * 3);
    memcpy(jac->data + 6, df3, sizeof(double) * 3);
    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = 3;
    jac->ncol = 3;
}


void impf_subrt_jacob_mfn(
    const impf_t_mfn_array fs, const unsigned int m, 
    const double * x, const unsigned int n, 
    double * buf, impf_t_matrix * jac)
{
    assert(NULL != fs);
    assert(NULL != x);
    assert(NULL != buf);
    assert(NULL != jac);

    double h, two_h, twelf_h;
    double *point, *ret1, *ret2, *ret3, *ret4;
    unsigned int i, j;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    two_h = 2 * h;
    twelf_h = 12 * h;

    point = buf;
    ret1 = buf + n;
    ret2 = ret1 + m;
    ret3 = ret2 + m;
    ret4 = ret3 + m;

    memcpy(point, x, sizeof(double) * n); /* Copy x into point */

    for(j = 0; j<n; j++)  /* Column j */
    {
        *(point + j) -= two_h;
        fs(m, n, point, ret1);  /* f(x_j - 2h) */
        *(point + j) += h;
        fs(m, n, point, ret2);  /* f(x_j - h)  */
        *(point + j) += two_h;
        fs(m, n, point, ret3);  /* f(x_j + h)  */
        *(point + j) += h;
        fs(m, n, point, ret4);  /* f(x_j + 2h) */
        *(point + j) -= two_h;

        for(i = 0; i<m; i++)  /* Row i */
        {
            *(jac->data + i * n + j) = \
            (*(ret1 + i) - *(ret4 + i) + 8 * (*(ret3 + i) - *(ret2 + i))) / twelf_h;
        }
    }

    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = m;
    jac->ncol = n;
}


void impf_subrt_jacob(
    const impf_t_mfn_vector fs, const unsigned int m, 
    const impf_t_vector * x, 
    double * buf, impf_t_matrix * jac)
{
    assert(NULL != fs);
    assert(NULL != x);
    assert(NULL != buf);
    assert(NULL != jac);

    impf_t_vector ptvec, ret1_vec, ret2_vec, ret3_vec, ret4_vec;
    double h, two_h, twelf_h;
    double *point, *ret1, *ret2, *ret3, *ret4;
    unsigned int i, j, n;

    /* Numerical Settings: step length */
    h = impf_SQRT_DBL_EPSILON;

    /* 5-Point Formula */
    n = x->dim;
    two_h = 2 * h;
    twelf_h = 12 * h;

    point = buf;
    ptvec.data = point; ptvec.dim = n;

    ret1 = buf + n;
    ret2 = ret1 + m;
    ret3 = ret2 + m;
    ret4 = ret3 + m;
    ret1_vec.data = ret1; ret1_vec.dim = m;
    ret2_vec.data = ret2; ret2_vec.dim = m;
    ret3_vec.data = ret3; ret3_vec.dim = m;
    ret4_vec.data = ret4; ret4_vec.dim = m;

    memcpy(point, x->data, sizeof(double) * n); /* Copy x into point */

    for(j = 0; j<n; j++)  /* Column j */
    {
        *(point + j) -= two_h;
        fs(&ptvec, &ret1_vec);  /* f(x_j - 2h) */
        *(point + j) += h;
        fs(&ptvec, &ret2_vec);  /* f(x_j - h)  */
        *(point + j) += two_h;
        fs(&ptvec, &ret3_vec);  /* f(x_j + h)  */
        *(point + j) += h;
        fs(&ptvec, &ret4_vec);  /* f(x_j + 2h) */
        *(point + j) -= two_h;

        for(i = 0; i<m; i++)  /* Row i */
        {
            *(jac->data + i * n + j) = \
            (*(ret1 + i) - *(ret4 + i) + 8 * (*(ret3 + i) - *(ret2 + i))) / twelf_h;
        }
    }

    jac->major = IMPF_MAT_ROW_MAJOR;
    jac->nrow = m;
    jac->ncol = n;
}

