/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <assert.h>
#include <stddef.h>
#include <impf/impf_basic.h>


/**
 * \biref BLAS subroutine DDOT
 *
 * Inner product: X'Y.
 */
extern double
ddot_(
    const int *N,     /**< [in] the number of selected elements of vector X and Y. */
    const double *DX, /**< [in] the vector X. */
    const int *INCX,  /**< [in] the spacing of vector X. */
    const double *DY, /**< [in] the vector Y. */
    const int *INCY   /**< [in] the spacing of vector Y. */
    );

/**
 * \brief BLAS subroutine DGEMV
 *
 * y = alpha*A[^T]*x + beta*y.
 */
extern void
dgemv_(
    const char *TRANS,   /**< [in] transpose or not, 'T' or 'N'. */
    const int *M,        /**< [in] number of rows of A. */
    const int *N,        /**< [in] number of columns of A. */
    const double* ALPHA, /**< [in] the first scaling factor. */
    const double *A,     /**< [in] the matrix of dimension (LDA, N). */
    const int *LDA,      /**< [in] the leading dimension of A, >= max(1, M). */
    const double *X,     /**< [in] the vector X. */
    const int* INCX,     /**< [in] the spacing of X. */
    const double *BETA,  /**< [in] the second scaling factor. */
    double *Y,           /**< [in,out] the vector Y. */
    const int *INCY      /**< [in] the spacing of Y. */
    );

/**
 * \brief BLAS subroutine DGEMM
 *
 * C = alpha*A[^T]*B[^T] + beta*C: A is M by K matrix and B is K by N matrix.
 */
extern void
dgemm_(
    const char *TRANSA,  /**< [in] transpose A or not, 'T' or 'N'. */
    const char *TRANSB,  /**< [in] transpose B or not, 'T' or 'N'. */
    const int *M,        /**< [in] number of rows of A[^T] and C. */
    const int *N,        /**< [in] number of columns of B[^T]. */
    const int *K,        /**< [in] number of columns of A[^T] and the number of rows of B[^T]. */
    const double *ALPHA, /**< [in] the 1st scaling factor. */
    const double *A,     /**< [in] the matrix A. */
    const int *LDA,      /**< [in] the leading dimension of matrix A. */
    const double *B,     /**< [in] the matrix B. */
    const int *LDB,      /**< [in] the leading dimension of matrix B. */
    const double *BETA,  /**< [in] the 2nd scaling factor. */
    double *C,           /**< [in,out] the matrix C. */
    const int *LDC       /**< [in] the leading dimension of matrix C. */
    );


/*-----------------------------------------------------------------------------
 * Vector & Matrix Multiplications: Calling BLAS
 *---------------------------------------------------------------------------*/

void 
impf_subrt_vec_inner(
    const impf_t_vector * v1,
    const unsigned int spacing1,
    const impf_t_vector * v2,
    const unsigned int spacing2,
    double * re)
{
    assert(NULL != v1);
    assert(NULL != v2);

    int N1, N2, N, INCX, INCY;

    N1 = 1 + (v1->dim - 1) / spacing1;
    N2 = 1 + (v2->dim - 1) / spacing2;
    N = N1 <= N2 ? N1 : N2;
    INCX = spacing1;
    INCY = spacing2;
    *re = ddot_(&N, v1->data, &INCX, v2->data, &INCY);
}

void
impf_subrt_mv_mul(
    const double alpha,                /* [in] the 1st scaling factor. */
    const impf_t_matrix * mat,         /* [in] the matrix. */
    const impf_t_submatinfo * submat,  /* [in] sub matrix information. */
    const char trans,                  /* [in] transpose the [sub-]matrix or not, 'T' or 'N'. */
    const impf_t_vector * v1,          /* [in] the 1st vector. */
    const unsigned int spacing1,       /* [in] the spacing of the 1st vector. */
    const double beta,                 /* [in] the 2nd scaling factor. */
    impf_t_vector * v2,                /* [in, out] the 2nd vector. */
    const unsigned int spacing2)       /* [in] the spacing of the 2nd vector. */
{
    assert(NULL != mat);
    assert(NULL != v1);
    assert(NULL != v2);

    if(IMPF_MAT_COL_MAJOR == mat->major)
    {
        double * A;
        int LDA, INCX, INCY, M, N;
        
        INCX = spacing1;
        INCY = spacing2;
        LDA = mat->nrow;

        /* Subtract sub-matrix */
        impf_submat_subtract(mat, submat, trans, &M, &N, &A);

        /* Calling BLAS function */
        dgemv_(&trans, &M, &N, &alpha, A, &LDA, v1->data, &INCX, &beta, v2->data, &INCY);
    }
    else /* Row major */
    {
        impf_t_matrix A_trans = {
            mat->data, mat->ncol, mat->nrow, IMPF_MAT_COL_MAJOR
        };
        char transA_trans = 'T' == trans ? 'N' : 'T';
        
        if(NULL != submat)
        {
            impf_t_submatinfo submat_trans = {
                submat->jloc, submat->iloc, submat->ncol, submat->nrow
            };
            impf_subrt_mv_mul(alpha, &A_trans, &submat_trans, transA_trans, v1, spacing1, beta, v2, spacing2);
        }
        else
        {
            impf_subrt_mv_mul(alpha, &A_trans, NULL, transA_trans, v1, spacing1, beta, v2, spacing2);
        }
    }
}

void
impf_subrt_mm_mul(
    const double alpha,                /* [in] the 1st scaling factor. */
    const impf_t_matrix * A,           /* [in] the matrix A: A[^T] of dim M by K. */
    const impf_t_submatinfo * submatA, /* [in] the sub-matrix information of A. */
    const char transA,                 /* [in] transpose the [sub-]matrix A or not, 'T' or 'N'. */
    const impf_t_matrix * B,           /* [in] the matrix B: B[^T] of dim K by N. */
    const impf_t_submatinfo * submatB, /* [in] the sub-matrix information of B. */
    const char transB,                 /* [in] transpose the [sub-]matrix B or not, 'T' or 'N'. */
    const double beta,                 /* [in] the 2nd scaling factor. */
    impf_t_matrix * C,                 /* [in, out] the matrix C: C[^T] of dim M by N. */
    const impf_t_submatinfo * submatC, /* [in] the sub-matrix information of C. */
    double * buf)                      /* a double array of length M * N if C is row major and beta is nonzero. */
{
    assert(NULL != A);
    assert(NULL != B);
    assert(NULL != C);

    if(IMPF_MAT_COL_MAJOR == A->major
    && IMPF_MAT_COL_MAJOR == B->major
    && IMPF_MAT_COL_MAJOR == C->major)
    {
        double * matA, * matB, * matC;
        int M, K, N, MA, MC, KA, KB, NB, NC, LDA, LDB, LDC;

        LDA = A->nrow;
        LDB = B->nrow;
        LDC = C->nrow;

        /* Sub-matrix A[^T]: M by K */
        impf_submat_subtract(A, submatA, transA, &MA, &KA, &matA);

        /* Sub-matrix B[^T]: K by N */
        impf_submat_subtract(B, submatB, transB, &KB, &NB, &matB);

        /* Sub-matrix C: M by N */
        impf_submat_subtract(C, submatC, 'N', &MC, &NC, &matC);

        M = MA <= MC ? MA : MC;
        K = KA <= KB ? KA : KB;
        N = NB <= NC ? NB : NC;
        dgemm_(&transA, &transB, &M, &N, &K, &alpha, matA, &LDA, matB, &LDB, &beta, matC, &LDC);
    }

    if(IMPF_MAT_ROW_MAJOR == A->major)
    {
        impf_t_matrix A_trans = {
            A->data, A->ncol, A->nrow, IMPF_MAT_COL_MAJOR
        };
        char transA_trans = 'T' == transA ? 'N' : 'T';

        if(NULL != submatA)
        {
            impf_t_submatinfo submatA_trans = {
                submatA->jloc, submatA->iloc, submatA->ncol, submatA->nrow
            };
            impf_subrt_mm_mul(alpha, &A_trans, &submatA_trans, transA_trans, B, submatB, transB, beta, C, submatC, buf);
        }
        else
        {
            impf_subrt_mm_mul(alpha, &A_trans, NULL, transA_trans, B, submatB, transB, beta, C, submatC, buf);
        }
    }
    else if(IMPF_MAT_ROW_MAJOR == B->major)
    {
        impf_t_matrix B_trans = {
            B->data, B->ncol, B->nrow, IMPF_MAT_COL_MAJOR
        };
        char transB_trans = 'T' == transB ? 'N' : 'T';

        if(NULL != submatB)
        {
            impf_t_submatinfo submatB_trans = {
                submatB->jloc, submatB->iloc, submatB->ncol, submatB->nrow
            };
            impf_subrt_mm_mul(alpha, A, submatA, transA, &B_trans, &submatB_trans, transB_trans, beta, C, submatC, buf);
        }
        else
        {
            impf_subrt_mm_mul(alpha, A, submatA, transA, &B_trans, NULL, transB_trans, beta, C, submatC, buf);
        }
    }
    else if(IMPF_MAT_ROW_MAJOR == C->major)
    {
        if(beta == 0)
        {
            C->major = IMPF_MAT_COL_MAJOR;
        }
        else
        {
            impf_mat_transmajor(C, buf);
        }
        impf_subrt_mm_mul(alpha, A, submatA, transA, B, submatB, transB, beta, C, submatC, buf);
    }
}


