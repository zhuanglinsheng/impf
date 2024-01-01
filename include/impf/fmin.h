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

enum impf_VariableType {
	impf_REAL,
	impf_INTEGER,
	impf_BINARY
};

enum impf_BoundType {
	impf_FR,	/* free */
	impf_UP,	/* upper bounded */
	impf_LO,	/* lower bounded */
	impf_BS		/* bounded from both sides */
};

enum impf_ConstraintType {
	impf_EQ,
	impf_GE,
	impf_LE
};

struct impf_VariableBound {
	char name[16];
	double lb;
	double ub;
	enum impf_BoundType b_type;
	enum impf_VariableType v_type;
};

struct impf_LinearConstraint {
	char name[16];
	double * coef;
	double rhs;
	enum impf_ConstraintType type;
};

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
