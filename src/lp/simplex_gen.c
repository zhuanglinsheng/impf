/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/fmin_lp.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

static int stdlpf_alloc(const int M, const int N, double **obj2, double **x2, double **coef2,
			struct impf_LinearConstraint **constraints2)
{
	*obj2 = NULL;
	*x2 = NULL;
	*coef2 = NULL;
	*constraints2 = NULL;

	*obj2 = impf_malloc(N * sizeof(double));
	if (*obj2 == NULL)
		return impf_EXIT_FAILURE;
	*x2 = impf_malloc(N * sizeof(double));
	if (*x2 == NULL) {
		impf_free(*obj2);
		return impf_EXIT_FAILURE;
	}
	*coef2 = impf_malloc(M * N * sizeof(double));
	if (*coef2 == NULL) {
		impf_free(*obj2);
		impf_free(*x2);
		return impf_EXIT_FAILURE;
	}
	*constraints2 = impf_malloc(M * sizeof(struct impf_LinearConstraint));
	if (*constraints2 == NULL) {
		impf_free(*obj2);
		impf_free(*x2);
		impf_free(*coef2);
		return impf_EXIT_FAILURE;
	}
	return impf_EXIT_SUCCESS;
}

static void stdlpf_free(double *obj2, double *x2, double *coef2, struct impf_LinearConstraint *constraints2)
{
	if (obj2)
		impf_free(obj2);
	if (x2)
		impf_free(x2);
	if (coef2)
		impf_free(coef2);
	if (constraints2)
		impf_free(constraints2);
}

/* Get the size of standard form LP
 *
 * Rules:
 *	1. free variable x = y1 - y2
 *	2. "x <= ub" will be added to constraints
 */
static void stdlpf_size(const struct impf_VariableBound *bounds, const int m, const int n, int *_M, int *_N)
{
	int j;
	*_M = m;
	*_N = n;

	for (j = 0; j < n; j++) {
		const struct impf_VariableBound *bd = bounds + j;

		if (impf_BOUND_T_FR == bd->b_type)
			(*_N)++;
		if (impf_BOUND_T_UP == bd->b_type || impf_BOUND_T_BS == bd->b_type)
			(*_M)++;
	}
}

/* Variable transformation 0: "xj <= ub"
 */
static void lp_transf_0(const double *objective, const struct impf_LinearConstraint *constraints,
			const struct impf_VariableBound *bounds, const int m, const int _N, const int j,
			int *ctr_var, int *ctr_ubcons, double *obj2,
			double *coef2, struct impf_LinearConstraint *constraints2)
{
	int i, idx;
	const struct impf_VariableBound *bd;

	bd = bounds + j;
	obj2[*ctr_var] = objective[j];

	if (impf_BOUND_T_UP == bd->b_type || impf_BOUND_T_BS == bd->b_type) {
		idx = m + (*ctr_ubcons);
		constraints2[idx].coef = coef2 + idx * _N;
		constraints2[idx].rhs = bd->ub;
		constraints2[idx].type = impf_CONS_T_LE;
		constraints2[idx].coef[*ctr_var] = 1.;
		(*ctr_ubcons)++;
	}
	for (i = 0; i < m; i++)
		coef2[(*ctr_var) + i * _N] = (constraints + i)->coef[j];
	(*ctr_var)++;
}

/* Variable transformation 1: "free xj" => "xj = y1 - y2" with y1, y2 >=0
 */
static void lp_transf_1(const double *objective, const struct impf_LinearConstraint *constraints,
			const struct impf_VariableBound *bounds, const int m, const int _N, const int j,
			int *ctr_var, double *obj2, double *coef2)
{
	int i;

	if (impf_BOUND_T_FR != (bounds + j)->b_type)
		return;
	obj2[*ctr_var] = -objective[j];

	for (i = 0; i < m; i++)
		coef2[(*ctr_var) + i * _N] = -((constraints + i)->coef[j]);
	(*ctr_var)++;
}

/* Variable transformation 2: "xj >= lb" => "y >= 0"
 */
static void lp_transf_2(const double *objective, const struct impf_LinearConstraint *constraints,
			const struct impf_VariableBound *bounds, const int m, const int _M, const int j,
			double *obj2, double *obj_diff, struct impf_LinearConstraint *constraints2)
{
	int i;
	const struct impf_VariableBound *bd = bounds + j;

	if (impf_BOUND_T_LO != bd->b_type && impf_BOUND_T_BS != bd->b_type)
		return;
	*obj_diff += objective[j] * bd->lb;

	for (i = 0; i < _M; i++) {
		struct impf_LinearConstraint *cons = constraints2 + i;

		if (i < m)
			cons->rhs -= ((constraints + i)->coef[j]) * (bd->lb);
		else
			cons->rhs -= bd->lb;
	}
}

/* Transform original LP into standard form
 *
 * Note:
 *	Original LP: allow for more variable bounds
 *	Standard LP: x >= 0
 */
static void lp_transstd(const double *objective, const struct impf_LinearConstraint *constraints,
			const struct impf_VariableBound *bounds, const int m, const int n, const int _M, const int _N,
			double *obj2, double *obj_diff, double *coef2, struct impf_LinearConstraint *constraints2)
{
	int i, j;
	int ctr_var = 0;
	int ctr_ubcons = 0;

	for (i = 0; i < m; i++) {
		constraints2[i].coef = coef2 + i * _N;
		constraints2[i].rhs = (constraints + i)->rhs;
		constraints2[i].type = (constraints + i)->type;
	}
	for (j = 0; j < n; j++) {
		lp_transf_0(objective, constraints, bounds, m, _N, j, &ctr_var, &ctr_ubcons, obj2, coef2, constraints2);
		lp_transf_1(objective, constraints, bounds, m, _N, j, &ctr_var, obj2, coef2);
		lp_transf_2(objective, constraints, bounds, m, _M, j, obj2, obj_diff, constraints2);
	}
}

/* Recover original LP solution and value from the standard form
 */
static void retreive_ori_lp_sol(const struct impf_VariableBound *bounds, const int n,
				const double *x2, const double value2, const double obj_diff, double *x, double *value)
{
	int ctr_var;
	int j;

	ctr_var = 0;

	for (j = 0; j < n; j++) {
		const struct impf_VariableBound *bd = bounds + j;

		if (impf_BOUND_T_FR == bd->b_type) {
			x[j] = x2[ctr_var] - x2[ctr_var + 1];
			ctr_var++;
		} else if (impf_BOUND_T_LO == bd->b_type || impf_BOUND_T_BS == bd->b_type) {
			x[j] = x2[ctr_var] + bd->lb;
		} else {
			x[j] = x2[ctr_var];
		}
		ctr_var++;
	}
	*value = value2 + obj_diff;
}

int impf_lp_simplex(const double *objective, const struct impf_LinearConstraint *constraints,
		    const struct impf_VariableBound *bounds,
		    const int m, const int n, const char *criteria, const int niter,
		    double *x, double *value, int *code)
{
	int _M, _N;
	double *obj2, *x2, *coef2;
	struct impf_LinearConstraint *constraints2;
	double value2 = 0, obj_diff = 0; /* value = value2 + obj_diff */

	assert(objective != NULL);
	assert(constraints != NULL);
	assert(x != NULL);
	assert(value != NULL);
	assert(code != NULL);

	if (bounds == NULL)
		return impf_lp_simplex_std(objective, constraints, m, n, criteria, niter, x, value, code);
	stdlpf_size(bounds, m, n, &_M, &_N);

	if (stdlpf_alloc(_M, _N, &obj2, &x2, &coef2, &constraints2) == impf_EXIT_FAILURE) {
		*code = impf_MemoryAllocError;
		return impf_EXIT_FAILURE;
	}
	impf_memset(coef2, 0., _M * _N);
	lp_transstd(objective, constraints, bounds, m, n, _M, _N, obj2, &obj_diff, coef2, constraints2);
	if (impf_lp_simplex_std(obj2, constraints2, _M, _N, criteria, niter, x2, &value2, code) == impf_EXIT_SUCCESS) {
		retreive_ori_lp_sol(bounds, n, x2, value2, obj_diff, x, value);
		stdlpf_free(obj2, x2, coef2, constraints2);
		*code = impf_Success;
		return impf_EXIT_SUCCESS;
	}
	stdlpf_free(obj2, x2, coef2, constraints2);
	return impf_EXIT_FAILURE; /* error code already updated */
}
