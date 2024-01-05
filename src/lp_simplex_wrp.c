/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/lp.h>

int impf_lp_simplex_wrp(const struct impf_Model_LP *model, const char *criteria, const int niter,
			double *x, double *value, int *code)
{
	int m = model->m;
	int n = model->n;
	double *obj = model->objective;
	struct impf_LinearConstraint *cons = model->constraints;
	struct impf_VariableBound *bounds = model->bounds;

	return impf_lp_simplex(obj, cons, bounds, m, n, criteria, niter, x, value, code);
}
