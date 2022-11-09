/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/**************************************************************************//**
 * \file impf_basic_blas.h
 * \brief Definitions of commonly used data structures and related methods.
 *****************************************************************************/

#ifndef __IMPF_BASIC_BLAS_H__
#define __IMPF_BASIC_BLAS_H__

#include <stdio.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_blas_structs 
 *****************************************************************************/
/**@{*/

/**
 * \brief Real vector.
 *
 * Eample file `eg_dts_impf_t_vec.c`:
 * \include eg_dts_impf_t_vec.c
 * \include eg_dts_impf_t_vec.out
 */
typedef struct
{
    double * data;      /**< \brief Where the datas are stored */
    unsigned int dim;   /**< \brief The dimensional of the vector */
} impf_t_vector;

/**
 * \brief Matrix data storage order
 */
typedef enum
{
    IMPF_MAT_ROW_MAJOR, /**< \brief Matrix is stored by row major    */
    IMPF_MAT_COL_MAJOR, /**< \brief Matrix is stored by column major */
} impf_t_matrix_major;

/**
 * \brief Real matrix.
 *
 * Example file `eg_dts_impf_t_mat.c`
 * \include eg_dts_impf_t_mat.c
 * \include eg_dts_impf_t_mat.out
 */
typedef struct
{
    double * data;              /**< \brief Where the data are stored. */
    unsigned int nrow;          /**< \brief Row number. */
    unsigned int ncol;          /**< \brief Column number. */
    impf_t_matrix_major major;  /**< \brief Data stored by row (1) or by column (0) */
} impf_t_matrix;

/**
 * \brief Submatrix Information.
 */
typedef struct
{
    unsigned int iloc;  /**< \brief the starting row index. */
    unsigned int jloc;  /**< \brief the starting column index. */
    unsigned int nrow;  /**< \brief the number of rows. */
    unsigned int ncol;  /**< \brief the number of columns. */
} impf_t_submatinfo;

/**
 * \brief Set all elements of matrix to be zeros.
 *
 * \details
 * Example: file `eg_dts_impf_mat_zeros.c`
 * \include eg_dts_impf_mat_zeros.c
 * \include eg_dts_impf_mat_zeros.out
 */
void
impf_mat_zeros(impf_t_matrix * mat);

/**
 * \brief Set all elements of matrix to be ones.
 *
 * \details
 * Example: file `eg_dts_impf_mat_ones.c`
 * \include eg_dts_impf_mat_ones.c
 * \include eg_dts_impf_mat_ones.out
 */
void
impf_mat_ones(impf_t_matrix * mat);

/**
 * \brief Set the principle elements of matrix to be one.
 *
 * \details
 * Example: file `eg_dts_impf_mat_eyes.c`
 * \include eg_dts_impf_mat_eyes.c
 * \include eg_dts_impf_mat_eyes.out
 */
void
impf_mat_eyes(impf_t_matrix * mat);

/**
 * \brief Transpose a matrix.
 *
 * \details Example: file `eg_dts_impf_mat_transpose.c`
 * \include eg_dts_impf_mat_transpose.c
 * \include eg_dts_impf_mat_transpose.out
 */
void
impf_mat_transpose(impf_t_matrix * mat);

/**
 * \brief Change the major of a matrix.
 * \param [in,out] mat the matrix to be transmajored.
 * \param buf a double array of the length mat->ncol * mat->nrow.
 *
 * \details Example: file `eg_dts_impf_mat_transmajor.c`
 * \include eg_dts_impf_mat_transmajor.c
 * \include eg_dts_impf_mat_transmajor.out
 */
void
impf_mat_transmajor(impf_t_matrix * mat, double * buf);

/**
 * \brief Change the major of a matrix (just for testing).
 * \param [in,out] mat the matrix to be transmajored.
 * \param buf a char array of the length mat->ncol * mat->nrow.
 *
 * \details In-place matrix transpose \cite cate1977algorithm. Generally, the
 * inclass matrix transpose is slower than using auxiliary memory.
 *
 * Also see: https://en.wikipedia.org/wiki/In-place_matrix_transposition
 *
 * Example: file `eg_dts_impf_mat_transmajor_inplace.c`
 * \include eg_dts_impf_mat_transmajor_inplace.c
 * \include eg_dts_impf_mat_transmajor_inplace.out
 */
void
impf_mat_transmajor_inplace(impf_t_matrix * mat, char * buf);

/**
 * \brief Subtract the submatrix from big matrix.
 *
 * \param [in] mat    the outer matrix.
 * \param [in] submat the submatrix info.
 * \param [in] trans  transpose the submatrix or not, 'T' or 'N'.
 * \param [out] M      number of rows in submat[^T].
 * \param [out] N      number of columns of submat[^T].
 * \param [out] data   starting point of submat.
 *
 * Example file `eg_dts_impf_submat_subtract.c`
 * \include eg_dts_impf_submat_subtract.c
 * \include eg_dts_impf_submat_subtract.out
 */
void
impf_submat_subtract(
    const impf_t_matrix * mat,
    const impf_t_submatinfo * submat, const char trans,
    int * M, int * N, double ** data);

/**
 * \brief Print a impf_t_vector in scientific notation.
 * \param [in] vec the impf_t_vector to be displayed.
 */
void
impf_print_vec_sci(const impf_t_vector * vec);

/**
 * \brief Print a impf_t_vector in scientific notation to file `out`.
 * \param [in,out] out the file pointer of the output
 * \param [in] vec the impf_t_vector to be displayed.
 */
void
impf_fprint_vec_sci(FILE * out, const impf_t_vector * vec);

/**
 * \brief Print a impf_t_matrix in scientific notation.
 * \param [in] mat the impf_t_matrix to be displayed.
 * \param [in] width the with of the printing panel, >= 40.
 * \param [in] foot display the size of \p mat. 0 or 1.
 */
void
impf_print_mat_sci(
    const impf_t_matrix * mat, const unsigned int width, const unsigned foot);

/**
 * \brief Print a impf_t_matrix in scientific notation to file `out`.
 * \param [in,out] out the file pointer of the output
 * \param [in] mat the impf_t_matrix to be displayed.
 * \param [in] width the with of the printing panel, >= 40.
 * \param [in] foot display the size of \p mat. 0 or 1.
 */
void
impf_fprint_mat_sci(
    FILE * out,
    const impf_t_matrix * mat, const unsigned int width, const unsigned foot);
/**@}*/

/**
 * \addtogroup module_blas_distances
 * \details
 * Norms and distances are implemented by calling BLAS.
 * @{
 */
/**************************************************************************//**
 * \ingroup module_blas_distances
 * \name Norms
 * Norm means the "size" of an impf_t_vector. Functions of L1-norm, L2-norm
 * and Linf-norm are provided. The definitions are
 *  \f[\begin{aligned}
 *      \|x\|_1 \equiv \sum_{i=1}^n|x_i|,\quad
 *      \|x\|_2 \equiv \left(\sum_{i=1}^nx_i^2\right)^\frac{1}{2},\quad
 *      \|x\|_\infty \equiv \max_{i\in\{1,2,...,n\}}|x_i|
 *  \end{aligned}\f]
 *
 * An example is given in `eg_dts_impf_df_norms.c`:
 * \include eg_dts_impf_df_norms.c
 * \include eg_dts_impf_df_norms.out
 *****************************************************************************/
/**@{*/

/** 
 * \brief 1-norm of a impf_t_vector.
 * \return a double float value, 1-norm of `vec`.
 * \param [in] vec the vector.
 * \param [in] spacing the spacing.
 */
double
impf_df_norm1(const impf_t_vector * vec, const unsigned int spacing);

/** 
 * \brief 2-norm of a impf_t_vector.
 * \return a double float value, 2-norm of `vec`.
 * \param [in] vec the vector.
 * \param [in] spacing the spacing.
 */
double
impf_df_norm2(const impf_t_vector * vec, const unsigned int spacing);

/** 
 * \brief Infinite-norm of a impf_t_vector.
 * \return a double float value, infinity norm of `vec`.
 * \param [in] vec the vector.
 * \param [in] spacing the spacing.
 */
double
impf_df_norminf(const impf_t_vector * vec, const unsigned int spacing);

/**
 * \brief The absolute value of a double float value.
 */
double
impf_df_abs(double x);

/**
 * \brief Infinite-norm of a 2d double array.
 */
double
impf_df_norminf_2d(const double * arr);

/**
 * \brief Infinite-norm of a 3d double array.
 */
double
impf_df_norminf_3d(const double * arr);
/**@}*/

/**************************************************************************//**
 * \ingroup module_blas_distances
 * \name Distances
 * The distance of two vectors are induced by norms: For any two vectors
 * \f$a, b\f$, their distance induced by Lp-norm is \f$\|a-b\|_p\f$. Concretly,
 *  \f[\begin{aligned}
 *      \|x-y\|_1 \equiv \sum_{i=1}^n|x_i-y_i|,\quad
 *      \|x-y\|_2 \equiv \left(\sum_{i=1}^n(x_i-y_i)^2\right)^\frac{1}{2},\quad
 *      \|x-y\|_\infty \equiv \max_{i\in\{1,2,...,n\}}|x_i-y_i|
 *  \end{aligned}\f]
 *
 * An example is given in `eg_dts_impf_df_distances.c`:
 * \include eg_dts_impf_df_distances.c
 * \include eg_dts_impf_df_distances.out
 *****************************************************************************/
/**@{*/

/**
 * \brief Absolute subtraction of two real values a and b, i.e., \f$|a-b|\f$.
 * \return a double float value, the norm absolute difference of `a` and `b`.
 * \param a [in] the 1st value.
 * \param b [in] the 2nd value.
 */
double
impf_df_abs_sub(const double a, const double b);

/**
 * \brief Norm-infinity distance of two 2D double float arrays.
 * \return a double float value, the norm-inf distance of `a` and `b`.
 * \param [in] a a nonempty pointer of a 2D double float array.
 * \param [in] b a nonempty pointer of a 2D double float array.
 */
double
impf_df_distance_norminf_2d(const double * a, const double * b);

/**
 * \brief Norm-infinity distance of two 3D double float arrays.
 * \return a double float value, the norm-inf distance of `a` and `b`.
 * \param [in] a a nonempty pointer of a 3D double float array.
 * \param [in] b a nonempty pointer of a 3D double float array.
 */
double
impf_df_distance_norminf_3d(const double * a, const double * b);

/**
 * \brief Norm-1 distance of two impf_t_vector. 
 * \details The dimension of the two impf_t_vector should be the same.
 * \return a double float value, the norm-1 distance of `a` and `b`.
 * \param [in] a the 1st vector.
 * \param [in] b the 2nd vector.
 * \param [in] spacing1 the spacing of the 1st vector.
 * \param [in] spacing2 the spacing of the 2nd vector.
 * \param [in] maxnums the maximum numbers of elements used in \p a and \p b.
 * \param [in] buf the calculation buffer, length = \p maxnums.
 */
double
impf_df_distance_norm1(
    const impf_t_vector * a, const impf_t_vector * b,
    const unsigned int spacing1, const unsigned int spacing2,
    const unsigned int maxnums, double * buf);

/**
 * \brief Norm-2 distance of two impf_t_vector. 
 * \details The dimension of the two impf_t_vector should be the same.
 * \return a double float value, the norm-2 distance of `a` and `b`.
 * \param [in] a the 1st vector.
 * \param [in] b the 2nd vector.
 * \param [in] spacing1 the spacing of the 1st vector.
 * \param [in] spacing2 the spacing of the 2nd vector.
 * \param [in] maxnums the maximum numbers of elements used in \p a and \p b.
 * \param [in] buf the calculation buffer, length = \p maxnums.
 */
double
impf_df_distance_norm2(
    const impf_t_vector * a, const impf_t_vector * b,
    const unsigned int spacing1, const unsigned int spacing2,
    const unsigned int maxnums, double * buf);

/**
 * \brief Norm-infinity distance of two impf_t_vector. 
 * \details The dimension of the two impf_t_vector should be the same.
 * \return a double float value, the norm-inf distance of `a` and `b`.
 * \param [in] a the 1st vector.
 * \param [in] b the 2nd vector.
 * \param [in] spacing1 the spacing of the 1st vector.
 * \param [in] spacing2 the spacing of the 2nd vector.
 * \param [in] maxnums the maximum numbers of elements used in \p a and \p b.
 * \param [in] buf the calculation buffer, length = \p maxnums.
 */
double
impf_df_distance_norminf(
    const impf_t_vector * a, const impf_t_vector * b,
    const unsigned int spacing1, const unsigned int spacing2,
    const unsigned int maxnums, double * buf);
/**@}*/
/** @} */

/**************************************************************************//**
 * \addtogroup module_blas_operations 
 * \details
 *
 * There are three basic (row) operations for matrices:
 *
 * 1. Swapping two rows;
 * 2. Scaling a row by a scalar;
 * 3. Plusing on a row another scaled row.
 *
 * Methods of vector based and matrix based versions of basic (row) operations
 * are defined below:
 *****************************************************************************/
/**@{*/

/**
 * \brief Swap two vectors.
 * \param [in,out] v1 the 1st vector whose selected elements are to be swaped.
 * \param [in,out] v2 the 2nd vector whose selected elements are to be swaped.
 * \param [in] spacing1 the spacing of \p v1.
 * \param [in] spacing2 the spacing of \p v2.
 */
void
impf_subrt_vec_swap(
    impf_t_vector * v1, impf_t_vector * v2,
    const unsigned int spacing1, const unsigned int spacing2);

/**
 * \brief Swap two rows of a matrix.
 * \param [in,out] mat the matrix whose selected elements in i-th and j-th rows are to be swapped.
 * \param [in] i the i-th row.
 * \param [in] j the j-th row.
 * \param [in] ispacing the spacing of i-th row.
 * \param [in] jspacing the spacing of j-th row.
 */
void
impf_subrt_mat_rowswap(
    impf_t_matrix * mat,
    const unsigned int i, const unsigned int j,
    const unsigned int ispacing, const unsigned int jspacing);

/**
 * \brief Vector multiplies a scalar: re = a * vec.
 * \param [in] vec the vector whose selected elements are to be scaled.
 * \param [in] a the scaling factor.
 * \param [in] spacing the spacing.
 * \param [out] re the scaled vector.
 */
void
impf_subrt_vec_scale(
    const impf_t_vector * vec, const double a,
    const unsigned int spacing, impf_t_vector * re);

/**
 * \brief Vector multiplied by a scalar: v = x * v.
 * \param [in,out] vec the vector whose selected elements are to be scaled.
 * \param [in] a the scaling factor.
 * \param [in] spacing the spacing.
 */
void
impf_subrt_vec_selfscale(
    impf_t_vector * vec, const double a, const unsigned int spacing);

/**
 * \brief Matrix whose i-th row multiplied by a scalar: Mi = Mi * a.
 * \param [in,out] mat the matrix whose selected elements in i-th row are to be scaled by \p a.
 * \param [in] a the scaling factor.
 * \param [in] i the i-th row.
 * \param [in] spacing the spacing.
 */
void impf_subrt_mat_rowscale(
    impf_t_matrix * mat,
    const double a, const unsigned int i, const unsigned int spacing);

/**
 * \brief Vector adds a vector multiplied by a scalar: re = a * v1 + v2.
 * \param [in] v1 the vector to be added by a scaled vector on selected elements.
 * \param [in] a the scaling factor.
 * \param [in] v2 the other vector to be added.
 * \param [in] spacing1 the spacing of \p v1.
 * \param [in] spacing2 the spacing of \p v2.
 * \param [out] re the the result.
 */
void
impf_subrt_vec_plusscaled(
    const impf_t_vector * v1, const double a, const impf_t_vector * v2,
    const unsigned int spacing1, const unsigned int spacing2, impf_t_vector * re);

/**
 * \brief Vector added by a vector multiplied by a scalar: v2 = a * v1 + v2.
 * \param [in] v1 the vector to be added by a scaled vector on selected elements.
 * \param [in] a the scaling factor.
 * \param [in] spacing1 the spacing of \p v1.
 * \param [in] spacing2 the spacing of \p v2.
 * \param [in,out] v2 the other vector, to be added by a * v1.
 */
void impf_subrt_vec_selfplusscaled(
    const impf_t_vector * v1, const double a,
    const unsigned int spacing1, const unsigned int spacing2, impf_t_vector * v2);

/**
 * \brief The j-th row of a matrix is added by the scaled i-th row: Mj = a * Mi + Mj.
 * \param [in,out] mat the matrix whose j-th row will be added by the scaled i-th row on selected elements.
 * \param [in] a the scaling factor.
 * \param [in] i the i-th row.
 * \param [in] j the j-th row.
 * \param [in] ispacing the spacing of the \p i-th row.
 * \param [in] jspacing the spacing of the \p j-th row.
 */
void impf_subrt_mat_rowplusscaled(
    impf_t_matrix * mat, const double a,
    const unsigned int i, const unsigned int j,
    const unsigned int ispacing, const unsigned int jspacing);
/**@}*/


/**************************************************************************//**
 * \addtogroup module_blas_multiply 
 *
 * \details 
 * 
 * Methods of inner product of vectors, matrix-vector product and 
 * matrix multiplications are defined here.
 *****************************************************************************/
/**@{*/

/**
 * \brief Inner product of two vectors.
 * \param [in] v1 the 1st vector.
 * \param [in] spacing1 the spacing of \p v1.
 * \param [in] v2 the 2nd vector.
 * \param [in] spacing2 the spacing of \p v2.
 * \param [out] re the inner product of v1 and v2.
 */
void impf_subrt_vec_inner(
    const impf_t_vector * v1, const unsigned int spacing1,
    const impf_t_vector * v2, const unsigned int spacing2, double * re);

/**
 * \brief Matrix-vector multiplication: v2 = alpha * M[^T] * v1 + beta * v2.
 * \param [in] alpha the 1st scaling factor.
 * \param [in] mat the matrix.
 * \param [in] submat to replace the mat with a sub-matrix of mat.
 * \param [in] trans transpose the [sub-]matrix or not, 'T' or 'N'.
 * \param [in] v1 the 1st vector.
 * \param [in] spacing1 the spacing of the 1st vector.
 * \param [in] beta the 2nd scaling factor.
 * \param [in,out] v2 the 2nd vector.
 * \param [in] spacing2 the spacing of the 2nd vector.
 */
void
impf_subrt_mv_mul(
    const double alpha,
    const impf_t_matrix * mat, const impf_t_submatinfo * submat, const char trans,
    const impf_t_vector * v1, const unsigned int spacing1,
    const double beta,
    impf_t_vector * v2, const unsigned int spacing2);

/**
 * \brief Matrix-matrix multiplication: C = alpha * A[^T] * B[^T] + beta * C.
 * \param [in] alpha the 1st scaling factor.
 * \param [in] A the matrix A: A[^T] of dim M by K.
 * \param [in] submatA to replace A with a sub-matrix of A.
 * \param [in] transA transpose the [sub-]matrix A or not, 'T' or 'N'.
 * \param [in] B the matrix B: B[^T] of dim K by N.
 * \param [in] submatB to replace B with a sub-matrix of B.
 * \param [in] transB transpose the [sub-]matrix B or not, 'T' or 'N'.
 * \param [in] beta the 2nd scaling factor.
 * \param [in,out] C the matrix C: C[^T] of dim M by N.
 * \param [in] submatC to replace C with a sub-matrix of C.
 * \param buf a double array of length M * N if C is row major and beta is nonzero.
 */
void
impf_subrt_mm_mul(
    const double alpha,
    const impf_t_matrix * A, const impf_t_submatinfo * submatA, const char transA,
    const impf_t_matrix * B, const impf_t_submatinfo * submatB, const char transB,
    const double beta,
    impf_t_matrix * C, const impf_t_submatinfo * submatC, double * buf);
/**@}*/

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_BASIC_BLAS_H__ */
