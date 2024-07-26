/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef __IMPF__MAGICS_H__
#define __IMPF__MAGICS_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

/*******************************************************************************
 * Magic numbers used in linear algebra and numerical differentiation
 ******************************************************************************/

#define __impf_IDF_LINALG_MAJ_ZERO__		1e-14
#define __impf_IDF_LINALG_DET_ZERO__		1e-9

/* The step of first-order differentiation in five-point formula
 * See "First order numerical differentiation family" in <diff.h>
 */
#define __impf_DIFF_H__				6e-4

/*******************************************************************************
 * Magic numbers uses in simplex algorithm
 *
 * Note: The following magics can be grouped to three types:
 * 	0. IDF, identifier
 *	1. CTR, controller
 *	2. CHC, checker
 * Identifier will not be directly used in algorithm
 * Controller affects the algorithm iterations
 * Checker only do some validity checks, not affecting the iterations
 ******************************************************************************/

/* Identifier of non-zero beta */
#define __impf_IDF_SPLX_ZEROS_BETA__		1e-9

/* Checker of the general checking "LP is optimal" */
#define __impf_CTR_SPLX_OPTIMAL__		1e-9

/* Checker of the checking "LP is feasible" */
#define __impf_CHC_SPLX_FEASIBLE__		1e-5

/* Checker of the checking "LP is degenerated" */
#define __impf_CHC_SPLX_DEGENERATED__		1e-12

/* Controller for pivot leaving rule */
#define __impf_CTR_SPLX_PIV_LEV__		1e-15

/* Controllers for Bland's rule */
#define __impf_CTR_SPLX_BLAND_EPS__		1e-6
#define __impf_CTR_SPLX_BLAND_EPS_MIN__		__impf_IDF_SPLX_ZEROS_BETA__

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
