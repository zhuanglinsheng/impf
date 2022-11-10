/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/*!****************************************************************************
 * \file impf_root.h
 * \brief Root Finding Algorithms
 *****************************************************************************/

#ifndef __IMPF_ROOT_H__
#define __IMPF_ROOT_H__

#include <impf/impf_basic.h>
#include <impf/impf_ndiff.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/**************************************************************************//**
 * \addtogroup module_root
 * \details 
 * To be added...
 * 
 * @{
 *****************************************************************************/

/**************************************************************************//**
 * \name Gauss-Newton Algorithms
 * \details This method search the unconstraint root of a function 
 * \f$f:\mathbb{R}^m\to\mathbb{R}^n\f$. 
 * 
 * When \f$m = n\f$, the algorithm is degerated to a classical newton's method.
 * We update the guesses \f$x\f$ by 
 * \f[
 * x_{t+1} = x_t - J^{-1}(x_t)f(x_t)
 * \f]
 * 
 * When \f$m > n\f$, it is natural to use \f$J^+ = (J^TJ)^{-1}J^T\f$ as a 
 * generalized Jacobian, and find a root in the nonlinear-least square sense. 
 * 
 * When \f$m < n\f$, 
 * 
 * APIs:
 *****************************************************************************/
/**@{*/

/**
 * \brief Finding the root of a 1 variant \ref realfunc "real function".
 * \param [in] f the real function of 1 variant.
 * \param [in,out] x the initial guess, will be updated and converge to the root.
 * \param [in] tol the precision of root.
 * \param [in] tolf the error of function value at root.
 * \param [in] niter the maximum iteration epochs.
 *
 * Example file `eg_root_newton_1f1.c`:
 * \include eg_root_newton_1f1.c
 * \include eg_root_newton_1f1.out
 */
impf_t_subrtstate
impf_subrt_root_newton_1f1(
    impf_t_f1 f, double * x,
    const double tol, const double tolf, const unsigned int niter);

/**
 * \brief Finding the root of two 2 variants \ref realfunc "real functions".
 * \param [in] f1 the real function of 2 variants.
 * \param [in] f2 the real function of 2 variants.
 * \param [in,out] x the initial guess, will be updated and converge to the root.
 * \param [in,out] y the initial guess, will be updated and converge to the root.
 * \param [in] tol the precision of root.
 * \param [in] tolf the error of function value at root.
 * \param [in] niter the maximum iteration epochs.
 *
 * Example file `eg_root_newton_2f2.c`:
 * \include eg_root_newton_2f2.c
 * \include eg_root_newton_2f2.out
 */
impf_t_subrtstate
impf_subrt_root_newton_2f2(
    impf_t_f2 f1, impf_t_f2 f2, double * x, double * y,
    const double tol, const double tolf, const unsigned int niter);

/**
 * \brief Finding the root of three 3 variants \ref realfunc "real functions".
 * \param [in] f1 the real function of 3 variants.
 * \param [in] f2 the real function of 3 variants.
 * \param [in] f3 the real function of 3 variants.
 * \param [in,out] x the initial guess, will be updated and converge to the root.
 * \param [in,out] y the initial guess, will be updated and converge to the root.
 * \param [in,out] z the initial guess, will be updated and converge to the root.
 * \param [in] tol the precision of root.
 * \param [in] tolf the error of function value at root.
 * \param [in] niter the maximum iteration epochs.
 *
 * Example file `eg_root_newton_3f3.c`:
 * \include eg_root_newton_3f3.c
 * \include eg_root_newton_3f3.out
 */
impf_t_subrtstate
impf_subrt_root_newton_3f3(
    impf_t_f3 f1, impf_t_f3 f2, impf_t_f3 f3,
    double * x, double * y, double * z,
    const double tol, const double tolf, const unsigned int niter);

/**@}*/

/** @} */ /*--- end of module_root ---*/

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_ROOT_H__ */
