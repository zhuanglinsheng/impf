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

#define __impf_IDF_LINALG_MAJ_ZERO__    1e-14
#define __impf_IDF_LINALG_DET_ZERO__    1e-9

/* The step of first-order differentiation in five-point formula
 * See "First order numerical differentiation family" in <diff.h>
 */
#define __impf_DIFF_H__                 6e-4

/*******************************************************************************
 * Magic numbers uses in simplex algorithm
 ******************************************************************************/

/* All `beta`-s less than this value is taken as zero in simplex tabular */
#define __impf_SPLX_ZEROS_BETA__        1e-12

/* Identifier of reasonable negative zero to be rounded */
#define __impf_SPLX_RHS_NZERO__         1e-9

/* Identifier of the general checking "LP is optimal" */
#define __impf_SPLX_OPTIMAL__           1e-9

/* Identifier of the checking "LP is feasible"  */
#define __impf_SPLX_FEASIBLE__          1e-5

/* Identifier of the checking "LP is degenerated" */
#define __impf_SPLX_DEGEN__             1e-10

/* Identifier of the checking "RHS of simplex table is nonnegative" */
#define __impf_SPLX_RHS_ZEROS__         __impf_SPLX_RHS_NZERO__

/* Identifier of the checking "the value of phase 1 simplex is nonnegative" */
#define __impf_SPLX_PHASE_1_NNVAL__     1e-5

/* Pivot leaving rule */
#define __impf_SPLX_PIVLEV_ZERO__       1e-9

/* Bland's rule */
#define __impf_SPLX_BLAND_EPS__         5e-4
#define __impf_SPLX_BLAND_EPS_MIN__     __impf_SPLX_ZEROS_BETA__

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
