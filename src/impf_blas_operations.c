/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <assert.h>
#include <string.h>
#include <impf/impf_basic_blas.h>


/**
 * \brief BLAS subroutine DSWAP.
 *
 * Swap vectors X and Y.
 */
void
dswap_(
    const int *N,    /**< [in] the length of selected elements in vector X. */
    double * DX,     /**< [in] the vector X. */
    const int *INCX, /**< [in] the spacing of vector X. */
    double *DY,      /**< [in,out] the vector Y. */
    const int *INCY  /**< [in] the spacing of vector Y. */
    );

/**
 * \brief BLAS subroutine DSCAL.
 *
 * X = a * X.
 */
void
dscal_(
    const int* N,     /**< [in] the length of selected elements in vector X. */
    const double* DA, /**< [in] the scalar a. */
    double* DX,       /**< [in,out] the vector X. */
    const int* INCX   /**< [in] the spacing of vector X. */
    );

/**
 * \brief BLAS subroutine DAXPY.
 *
 * Y = a * X + Y.
 */
void
daxpy_(
    const int* N,     /**< [in] the length of selected elements in vector X and Y. */
    const double* DA, /**< [in] the scalar a. */
    const double* DX, /**< [in] the vector X. */
    const int* INCX,  /**< [in] the spacing of vector X. */
    double* DY,       /**< [in] the vector Y. */
    const int* INCY   /**< [in] the spacing of vector Y. */
    );


/*-----------------------------------------------------------------------------
 * Vector & Matrix Operations: Calling BLAS
 *---------------------------------------------------------------------------*/

void 
impf_subrt_vec_swap(
    impf_t_vector * v1,         
    impf_t_vector * v2,         
    const unsigned int spacing1, 
    const unsigned int spacing2) 
{
    assert(NULL != v1);
    assert(NULL != v2);

    int N, N1, N2, inc1, inc2;
    N1 = 1 + (v1->dim - 1) / spacing1;
    N2 = 1 + (v2->dim - 1) / spacing2;
    N = N1 <= N2 ? N1 : N2;
    inc1 = spacing1;
    inc2 = spacing2;
    dswap_(&N, v1->data, &inc1, v2->data, &inc2);
}

void 
impf_subrt_mat_rowswap(
    impf_t_matrix * mat,      
    const unsigned int i,    
    const unsigned int j,     
    const unsigned int ispacing,
    const unsigned int jspacing)  
{
    assert(NULL != mat);

    if(IMPF_MAT_ROW_MAJOR == mat->major)
    {
        int N, N1, N2, inc_i, inc_j;

        inc_i = ispacing; 
        inc_j = jspacing; 
        N1 = 1 + (mat->ncol - 1) / inc_i;
        N2 = 1 + (mat->ncol - 1) / inc_j;
        N = N1 <= N2 ? N1 : N2;
        double * row_i = mat->data + i * mat->ncol;
        double * row_j = mat->data + j * mat->ncol;
        dswap_(&N, row_i, &inc_i, row_j, &inc_j);
    }
    if(IMPF_MAT_COL_MAJOR == mat->major)
    {
        double tmp;
        unsigned int iloc, jloc, idx_col_i = 0, idx_col_j = 0;
    SWAPING: 
        iloc = idx_col_i * (mat->nrow) + i;
        jloc = idx_col_j * (mat->nrow) + j;
        tmp = *(mat->data + iloc);
        *(mat->data + iloc) = *(mat->data + jloc);
        *(mat->data + jloc) = tmp;
        idx_col_i += ispacing;
        idx_col_j += jspacing;
        if(idx_col_i < mat->ncol && idx_col_j < mat->ncol) goto SWAPING;
    }
}

void 
impf_subrt_vec_scale(
    const impf_t_vector * vec,  
    const double a,             
    const unsigned int spacing, 
    impf_t_vector * re)          
{
    assert(NULL != vec);
    assert(NULL != re);

    int N, len, inc;
    len = vec->dim <= re->dim ? vec->dim : re->dim;
    N = 1 + (len - 1) / spacing;
    inc = spacing;
    memset(re->data, 0, sizeof(double) * re->dim);
    memcpy(re->data, vec->data, sizeof(double) * len);
    dscal_(&N, &a, re->data, &inc);
}

void 
impf_subrt_vec_selfscale(
    impf_t_vector * vec,        
    const double a,              
    const unsigned int spacing) 
{
    assert(NULL != vec);

    int N, inc;
    N = 1 + (vec->dim - 1) / spacing;
    inc = spacing;
    dscal_(&N, &a, vec->data, &inc);
}

void 
impf_subrt_mat_rowscale(
    impf_t_matrix * mat,         
    const double a,           
    const unsigned int i,        
    const unsigned int spacing) 
{
    assert(NULL != mat);

    if(IMPF_MAT_ROW_MAJOR == mat->major)
    {
        int N, inc;

        N = 1 + (mat->ncol - 1) / spacing;
        inc = spacing;
        dscal_(&N, &a, mat->data + i * N, &inc);
    }
    else
    {
        unsigned int idx_col = 0;
    SCALING:
        *(mat->data + idx_col * (mat->nrow) + i) *= a;
        idx_col += spacing;
        if(idx_col < mat->ncol) goto SCALING;
    }
}

void 
impf_subrt_vec_plusscaled(
    const impf_t_vector * v1,   
    const double a,              
    const impf_t_vector * v2,   
    const unsigned int spacing1, 
    const unsigned int spacing2, 
    impf_t_vector * re)          
{
    assert(NULL != v1);
    assert(NULL != v2);

    int incx, incy, N, N1, N2, len2;
    len2 = v2->dim <= re->dim ? v2->dim : re->dim;
    N1 = 1 + (v1->dim - 1) / spacing1;
    N2 = 1 + (len2 - 1) / spacing2;
    N = N1 <= N2 ? N1 : N2;
    incx = spacing1;
    incy = spacing2;
    memcpy(re->data, v2->data, sizeof(double) * len2);
    daxpy_(&N, &a, v1->data, &incx, re->data, &incy);
}

void 
impf_subrt_vec_selfplusscaled(
    const impf_t_vector * v1,     
    const double a,             
    const unsigned int spacing1, 
    const unsigned int spacing2,  
    impf_t_vector * v2)         
{
    assert(NULL != v1);
    assert(NULL != v2);

    int incx, incy, N, N1, N2;
    N1 = 1 + (v1->dim - 1) / spacing1;
    N2 = 1 + (v2->dim - 1) / spacing2;
    N = N1 <= N2 ? N1 : N2;
    incy = spacing2;
    incx = spacing1;
    daxpy_(&N, &a, v1->data, &incx, v2->data, &incy);
}

void impf_subrt_mat_rowplusscaled(
    impf_t_matrix * mat,          
    const double a,               
    const unsigned int i,         
    const unsigned int j,         
    const unsigned int ispacing,  
    const unsigned int jspacing)
{
    assert(NULL != mat);

    if(IMPF_MAT_ROW_MAJOR == mat->major)
    {
        int N, spacing, inc_i, inc_j;
        spacing = ispacing >= jspacing ? ispacing : jspacing;
        N = 1 + (mat->ncol - 1) / spacing;
        inc_i = ispacing;
        inc_j = jspacing;
        double * DX = mat->data + i * mat->ncol;  /* the i-th row. */
        double * DY = mat->data + j * mat->ncol;  /* the j-th row. */
        daxpy_(&N, &a, DX, &inc_i, DY, &inc_j);
    }
    else 
    {
        unsigned int iloc, jloc, idx_col_i = 0, idx_col_j = 0;
    PLUSSCALING:
        iloc = idx_col_i * mat->nrow + i;  /* row i, col idx_col */
        jloc = idx_col_j * mat->nrow + j;  /* row j, col idx_col */
        *(mat->data + jloc) += *(mat->data + iloc) * a;
        idx_col_i += ispacing;
        idx_col_j += jspacing;
        if(idx_col_i < mat->ncol && idx_col_j < mat->ncol) goto PLUSSCALING;
    }
}
