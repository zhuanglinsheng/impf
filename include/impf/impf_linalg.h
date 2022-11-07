/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**************************************************************************//**
 * \file impf_linalg.h
 * \brief A Minimum Implimentation of Linear Algebra Algorithms.
 *****************************************************************************/

#ifndef __IMPF_LAPACK_H__
#define __IMPF_LAPACK_H__

#include <impf/impf_basic.h>
#include <impf/impf_basic_blas.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_linalg
 * 
 * \details 
 * TO BE ADDED...
 * 
 * @{
 *****************************************************************************/

/**
 * \brief Solve the linear equation A * X = b where A is 2 by 2, b is 2 by 1.
 * \param [in] A11 the element A[1,1] of matrix A.
 * \param [in] A12 the element A[1,2] of matrix A.
 * \param [in] A21 the element A[2,1] of matrix A.
 * \param [in] A22 the element A[2,2] of matrix A.
 * \param [in] b1 the first element of vector B.
 * \param [in] b2 the second element of vector B.
 * \param [out] x the solution of Ax = b, a double array of the length 2.
 *
 * Example file `eg_linalg_lsolve_2A2_2b.c`:
 * \include eg_linalg_lsolve_2A2_2b.c
 * \include eg_linalg_lsolve_2A2_2b.out
 */
impf_t_subrtstate impf_subrt_lsolve_2A2_b2(
    const double A11, const double A12, const double A21, const double A22,
    const double b1, const double b2, double * x);

/**
 * \brief Solve the linear equation A * X = b where A is 3 by 3, b is 3 by 1.
 * \param [in] A11 the element A[1,1] of matrix A.
 * \param [in] A12 the element A[1,2] of matrix A.
 * \param [in] A13 the element A[1,3] of matrix A.
 * \param [in] A21 the element A[2,1] of matrix A.
 * \param [in] A22 the element A[2,2] of matrix A.
 * \param [in] A23 the element A[2,3] of matrix A.
 * \param [in] A31 the element A[3,1] of matrix A.
 * \param [in] A32 the element A[3,2] of matrix A.
 * \param [in] A33 the element A[3,3] of matrix A.
 * \param [in] b1 the 1st element of vector B.
 * \param [in] b2 the 2nd element of vector B.
 * \param [in] b3 the 3rd element of vector B.
 * \param [out] x the solution of Ax = b, a double array of the length 3.
 *
 * Example file `eg_linalg_lsolve_3A3_3b.c`:
 * \include eg_linalg_lsolve_3A3_3b.c
 * \include eg_linalg_lsolve_3A3_3b.out
 */
impf_t_subrtstate impf_subrt_lsolve_3A3_b3(
    const double A11, const double A12, const double A13,
    const double A21, const double A22, const double A23,
    const double A31, const double A32, const double A33,
    const double b1, const double b2, const double b3, double * x);

/**
 * \brief Solve the linear equation A[^T] * X = B[^T].
 */
impf_t_subrtstate impf_subrt_lsolve(
    const impf_t_matrix * A, const impf_t_submatinfo * submatA, const char transA,
    impf_t_matrix * X, int * IPIV,
    const impf_t_matrix * B, const impf_t_subrtstate * submatB, const char transB);



/** @}  */  /* end of module_linalg. */

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_LAPACK_H__ */
