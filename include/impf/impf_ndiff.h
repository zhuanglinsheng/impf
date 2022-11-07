/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/*!****************************************************************************
 * \file impf_ndiff.h
 * \brief Numerical Differentiation.
 *****************************************************************************/

#ifndef __IMPF_NDIFF_H__
#define __IMPF_NDIFF_H__

#include <impf/impf_basic.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_ndiff_grad
 *
 * \details 
 * The gradient vector of \f$f(x): \mathbb{R}^n\to\mathbb{R}\f$ is a 
 * \f$n\f$-dimensional double array. 
 *
 * The derivatives are calculated using a
 * [5-points method](https://en.wikipedia.org/wiki/Five-point_stencil):
 * Given a step length \f$h\f$, the derivative of real function
 * \f$f:\mathbb{R}\to\mathbb{R}\f$ at \f$x\f$ is calculated by
 * \f[
 *     f'(x) \approx \frac{f(x-2h)-f(x+2h)+8[f(x+h)-f(x-h)]}{h}
 * \f]
 * The error is of the degree \f$\mathcal{O}(h^5)\f$. The step length \f$h\f$
 * is chosen as \f$\sqrt{\epsilon}\f$, where \f$\epsilon\f$ is the double
 * epsilon defined in `<float.h>`.
 *****************************************************************************/
/**@{*/

/**
 * \brief Calculate the derivative of a 1-variant real function.
 * \param [in] f the differentiated function.
 * \param [in] x the evaluation point.
 * \param [out] diff the derivative of f at x.
 */
void impf_subrt_diff(const impf_t_f1 f, const double x, double * diff);

/**
 * \brief Calculate the gradient of a 2-variant real function.
 * \param [in] f the differentiated function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [out] grd the gradient of f at (x, y).
 * 
 * Example file `eg_diff_impf_grd_f2.c`
 * \include eg_diff_impf_grd_f2.c
 * ```
 * The gradient of f2 at (0.333333,0.714286) is (2.095238,2.095238)
 * ```
 */
void impf_subrt_gradient_f2(const impf_t_f2 f, const double x, const double y, double * grd);

/**
 * \brief Calculate the gradient of a 3-variant real function.
 * \param [in] f the differentiated function.
 * \param [in] x the 1st scalar of the evaluation point.
 * \param [in] y the 2nd scalar of the evaluation point.
 * \param [in] z the 3rd scalar of the evaluation point.
 * \param [out] grd the gradient of f at (x, y, z).
 */
void impf_subrt_gradient_f3(const impf_t_f3 f, const double x, const double y, const double z, double * grd);

/**
 * \brief Calculate the gradient of a n-variant real function.
 * \param [in] f the differentiated function.
 * \param [in] x the evaluation point.
 * \param [in] n the length of x.
 * \param buf the calclation buffer, a n-dimensional double array.
 * \param [out] grd the gradient of f at x.
 */
void impf_subrt_gradient_fn(const impf_t_fn f, const double * x, const unsigned int n, double * buf, double * grd);

/**
 * \brief Calculate the gradient of a n-variant real function.
 * \param [in] f the differentiated function.
 * \param [in] vec the evaluation point.
 * \param buf the calculation buffer, a n-dimensional double array.
 * \param [out] grd the gradient of f at vec.
 */
void impf_subrt_gradient(const impf_t_fn f, const impf_t_vector * vec, double * buf, double * grd);
/**@}*/


/**************************************************************************//**
 * \addtogroup module_ndiff_jacb
 *****************************************************************************/
/**@{*/

/**
 * \brief Calculate the Jacobian matrix of 2 functions of 1 variant.
 */
void impf_subrt_jacob_2f1(
    const impf_t_f1 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f1 f2,      /**< [in] the 2nd differentiated function. */ 
    const double x,          /**< [in] the evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$2\times1\f$, row major.  */
    );

/**
 * \brief Calculate the Jacobian matrix of 2 functions of 2 variants. 
 */
void impf_subrt_jacob_2f2(
    const impf_t_f2 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f2 f2,      /**< [in] the 2nd differentiated function. */ 
    const double x,          /**< [in] the 1st evaluation point. */
    const double y,          /**< [in] the 2nd evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$2\times2\f$, row major. */
    );

/**
 * \brief Calculate the Jacobian matrix of 2 functions of 3 variants. 
 */
void impf_subrt_jacob_2f3(
    const impf_t_f3 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f3 f2,      /**< [in] the 2nd differentiated function. */ 
    const double x,          /**< [in] the 1st evaluation point. */
    const double y,          /**< [in] the 2nd evaluation point. */
    const double z,          /**< [in] the 3rd evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$2\times3\f$, row major. */
    );

/**
 * \brief Calculate the Jacobian matrix of 3 functions of 1 variant.
 */
void impf_subrt_jacob_3f1(
    const impf_t_f1 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f1 f2,      /**< [in] the 2nd differentiated function. */ 
    const impf_t_f1 f3,      /**< [in] the 3rd differentiated function. */ 
    const double x,          /**< [in] the evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$3\times1\f$, row major.  */
    );

/**
 * \brief Calculate the Jacobian matrix of 3 functions of 2 variants. 
 */
void impf_subrt_jacob_3f2(
    const impf_t_f2 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f2 f2,      /**< [in] the 2nd differentiated function. */ 
    const impf_t_f2 f3,      /**< [in] the 3rd differentiated function. */ 
    const double x,          /**< [in] the 1st evaluation point. */
    const double y,          /**< [in] the 2nd evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$3\times2\f$, row major. */
    );

/**
 * \brief Calculate the Jacobian matrix of 3 functions of 3 variants. 
 */
void impf_subrt_jacob_3f3(
    const impf_t_f3 f1,      /**< [in] the 1st differentiated function. */ 
    const impf_t_f3 f2,      /**< [in] the 2nd differentiated function. */ 
    const impf_t_f3 f3,      /**< [in] the 3rd differentiated function. */ 
    const double x,          /**< [in] the 1st evaluation point. */
    const double y,          /**< [in] the 2nd evaluation point. */
    const double z,          /**< [in] the 3rd evaluation point. */
    impf_t_matrix * jac      /**< [out] the Jacobian matrix of dimension \f$3\times3\f$, row major. */
    );

/** 
 * \brief Calculate the Jacobian matrix of a impf_t_mfn_array function.
 * \param [in] mfn the differentiated function.
 * \param [in] m the dimension of return of \p mfn.
 * \param [in] x the evaluation point.
 * \param [in] n the length of x.
 * \param buf the calclation buffer, a (n+4m)-dimensional double array.
 * \param [out] jac the Jacobian matrix of dimension \f$m\times n\f$, row major.
 */
void impf_subrt_jacob_mfn(
    const impf_t_mfn_array mfn, const unsigned int m, const double * x, const unsigned int n, 
    double * buf, 
    impf_t_matrix * jac);

/** 
 * \brief Calculate the Jacobian matrix of a impf_t_mfn_vector function.
 * \param [in] mfn the differentiated function.
 * \param [in] m the dimension of return of \p mfn.
 * \param [in] x the evaluation point.
 * \param buf the calclation buffer, a (n+4m)-dimensional double array.
 * \param [out] jac the Jacobian matrix of dimension \f$m\times n\f$, row major.
 * \details 
 * The Jacobian matrix of \f$f(x): \mathbb{R}^n\to\mathbb{R}^m\f$ is a 
 * \f$m\times n\f$ matrix. 
 *
 * Example file `eg_diff_impf_jacob.c`. For a function
 * \f[
 * f_i(x_1, ..., x_n) = \sum^n_{j=1}(x_j-i-j)^2
 * \f]
 * The analytical Jacobian matrix is given by 
 * \f[
 * J_{i,j} = \frac{\partial f_i(x_1, ..., x_n)}{\partial x_j} = 2(x_j - i - j)
 * \f]
 * \include eg_diff_impf_jacob.c
 * \include eg_diff_impf_jacob.out
 */
void impf_subrt_jacob(
    const impf_t_mfn_vector mfn, const unsigned int m, const impf_t_vector * x,
    double * buf,              
    impf_t_matrix * jac);

/**@}*/

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_NDIFF_H__ */
