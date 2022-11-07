/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <assert.h>
#include <stddef.h>
#include <math.h>
#include <impf/impf_basic.h>


/**
 * \biref BLAS subroutine DASUM
 *
 * Sum of absolute values.
 */
extern double
dasum_(
    const int *N,     /**< [in] the number of selected elements of vector X. */
    const double *DX, /**< [in] the vector X. */
    const int *INCX   /**< [in] the spacing of vector X. */
    );

/**
 * \brief BLAS subroutine DNRM
 *
 * Norm 2.
 */
extern double
dnrm2_(const int *N, const double *DX, const int *INCX);

/**
 * \brief BLAS subroutine IDAMAX
 *
 * Location of absolute maximum value.
 */
extern int
idamax_(
    const int *N,     /**< [in] the number of selected elements of vector X. */
    const double *DX, /**< [in] the vector X. */
    const int *INCX   /**< [in] the spacing of vector X. */
    );

#ifndef __impf_MAX__
#define __impf_MAX__(a, b) ((a) > (b) ? (a) : (b))
#endif

/*-----------------------------------------------------------------------------
 * Norms: Calling BLAS.
 *---------------------------------------------------------------------------*/

double
impf_df_norm1(const impf_t_vector * vec, const unsigned int spacing)
{
    int incx, n;
    incx = spacing;
    n = 1 + (vec->dim - 1) / spacing;
    return dasum_(&n, vec->data, &incx);
}

double
impf_df_norm2(const impf_t_vector * vec, const unsigned int spacing)
{
    int incx, n;
    incx = spacing;
    n = 1 + (vec->dim - 1) / spacing;
    return dnrm2_(&n, vec->data, &incx);
}

double
impf_df_norminf(const impf_t_vector * vec, const unsigned int spacing)
{
    int incx, n, loc;
    incx = spacing;
    n = 1 + (vec->dim - 1) / spacing;
    loc = idamax_(&n, vec->data, &incx) - 1; /* Fortran is 1-indexed. */
    return fabs(*(vec->data + loc));
}

double
impf_df_abs(double x)
{
    return fabs(x);
}

double
impf_df_norminf_2d(const double * arr)
{
    double tmp1, tmp2;
    tmp1 = fabs(arr[0]);
    tmp2 = fabs(arr[1]);
    return __impf_MAX__(tmp1, tmp2);
}

double
impf_df_norminf_3d(const double * arr)
{
    double tmp1, tmp2, tmp3;
    tmp1 = fabs(arr[0]);
    tmp2 = fabs(arr[1]);
    tmp3 = fabs(arr[2]);
    return __impf_MAX__(__impf_MAX__(tmp1, tmp2), tmp3);
}

/*-----------------------------------------------------------------------------
 * Distances: 1d, 2d, 3d. Calling <math.h> function fabs
 *---------------------------------------------------------------------------*/

double 
impf_df_abs_sub(const double a, const double b)
{
    return fabs(a - b);
}

double 
impf_df_distance_norminf_2d(const double * a, const double * b)
{
    assert(NULL != a);
    assert(NULL != b);

    double tmp1, tmp2;

    tmp1 = fabs(*(a    ) - *(b    ));
    tmp2 = fabs(*(a + 1) - *(b + 1));
    return __impf_MAX__(tmp1, tmp2);
}

double 
impf_df_distance_norminf_3d(const double * a, const double * b)
{
    assert(NULL != a);
    assert(NULL != b);

    double tmp1, tmp2, tmp3;
    tmp1 = fabs(*(a    ) - *(b    ));
    tmp2 = fabs(*(a + 1) - *(b + 1));
    tmp3 = fabs(*(a + 2) - *(b + 2));
    return __impf_MAX__(__impf_MAX__(tmp1, tmp2), tmp3);
}

/*-----------------------------------------------------------------------------
 * Distances: Calling Norms.
 *---------------------------------------------------------------------------*/

double 
impf_df_distance_norm1(
    const impf_t_vector * a,    
    const impf_t_vector * b,      
    const unsigned int spacing1,  
    const unsigned int spacing2, 
    const unsigned int maxnums,  
    double * buf)                
{
    assert(NULL != a);
    assert(NULL != b);
    assert(NULL != buf);

    impf_t_vector b_minus_a = {buf, maxnums};
    impf_subrt_vec_plusscaled(a, -1, b, spacing1, spacing2, &b_minus_a);
    return impf_df_norm1(&b_minus_a, 1);
}

double 
impf_df_distance_norm2(
    const impf_t_vector * a,   
    const impf_t_vector * b,     
    const unsigned int spacing1, 
    const unsigned int spacing2, 
    const unsigned int maxnums,  
    double * buf)            
{
    assert(NULL != a);
    assert(NULL != b);
    assert(NULL != buf);

    impf_t_vector b_minus_a = {buf, maxnums};
    impf_subrt_vec_plusscaled(a, -1, b, spacing1, spacing2, &b_minus_a);
    return impf_df_norm2(&b_minus_a, 1);
}

double
impf_df_distance_norminf(
    const impf_t_vector * a, 
    const impf_t_vector * b,
    const unsigned int spacing1, 
    const unsigned int spacing2, 
    const unsigned int maxnums,  
    double * buf)               
{
    assert(NULL != a);
    assert(NULL != b);
    assert(NULL != buf);

    impf_t_vector b_minus_a = {buf, maxnums};
    impf_subrt_vec_plusscaled(a, -1, b, spacing1, spacing2, &b_minus_a);
    return impf_df_norminf(&b_minus_a, 1);
}
