/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

/*!****************************************************************************
 * \file impf_approx.h
 * \brief Using polynomials as bases to approximate the implicit function.
 *****************************************************************************/

#ifndef __IMPF_APPROX_H__
#define __IMPF_APPROX_H__

#include <impf/impf_basic.h>
#include <impf/impf_nls.h>

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */


/**************************************************************************//**
 * \addtogroup module_poly_approx
 *
 * @{
 *****************************************************************************/

typedef struct 
{
    unsigned int order;
    double * coefficients;
} impf_t_poly;



/**************************************************************************//**
 * \name Subroutines
 *****************************************************************************/
/**@{*/
 
/**@}*/




/**************************************************************************//**
 * \name Functions
 *****************************************************************************/
/**@{*/

/**
 * \brief Sampling polynomial bases on a cube.
 * \param m unsigned integer, number of implicit functions.
 * \param n unsigned integer, number of parameters of each impf.
 * \param poly_order unsigned integer, the order of polynomials.
 * \param lbs a nonempty double pointer, lower bounds of the cube.
 * \param ubs a nonempty double pointer, upper bounds of the cube.
 * \return the sampled matrix.
 */
    impf_t_matrix * impf_func_polybases_sampling_cube(
    const unsigned int m, 
    const unsigned int n, 
    const unsigned int poly_order, 
    const double * lbs, const double * ubs);
/**@}*/

/** @} */ /*--- end of module_poly_approx ---*/

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif /* __IMPF_APPROX_H__ */
