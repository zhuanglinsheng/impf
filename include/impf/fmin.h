/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include "utils.h"

#ifndef __IMPF_FMIN_H__
#define __IMPF_FMIN_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

#define impf_VAR_T_REAL	0
#define impf_VAR_T_INT	1
#define impf_VAR_T_BIN	2

#define impf_BOUND_T_FR	0	/* free */
#define impf_BOUND_T_UP	1	/* upper bounded */
#define impf_BOUND_T_LO	2	/* lower bounded */
#define impf_BOUND_T_BS	3	/* bounded from both sides */

#define impf_CONS_T_EQ	0
#define impf_CONS_T_GE	1
#define impf_CONS_T_LE	2

struct impf_VariableBound {
	char name[16];
	double lb;
	double ub;
	int b_type;
	int v_type;
};

struct impf_LinearConstraint {
	char name[16];
	double * coef;
	double rhs;
	int type;
};

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
