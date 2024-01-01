/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/lp.h>
#include <stdio.h>

/* LP Example
 * (from: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.linprog.html)
 *
 *         min      - x0 + 4 * x1
 *         s.t.  -3 * x0 +     x1 <=  6
 *                   -x0 - 2 * x1 >= -4
 *                             x1 >= -3
 *
 * The solution is (10, -3) and the optimal value is -22
 */
#define m 2        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-1., 4.};
double constraint_1_coef[] = {-3., 1.};
double constraint_2_coef[] = {-1., -2.};

struct impf_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  6., impf_LE },
	{ "", constraint_2_coef, -4., impf_GE },
};
struct impf_VariableBound bounds[] = {
	{ "x0", __impf_NINF__, __impf_INF__, impf_FR, impf_REAL },
	{ "x1",            -3, __impf_INF__, impf_LO, impf_REAL },
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	enum impf_ErrorCode code;
	int state = impf_lp_simplex(obj, constraints, bounds, m, n, "bland", 1000, x, &value, &code);

	printf("Error code = %u\n", code);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(value + 22) < 1e-8);
	assert(__impf_ABS__(x[0] - 10.) < 1e-8);
	assert(__impf_ABS__(x[1] + 3.) < 1e-8);
	return 0;
}
