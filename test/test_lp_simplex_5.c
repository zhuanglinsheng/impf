/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/fmin_lp.h>
#include <stdio.h>

/* LP Example
 * (from: https://www.mathworks.com/help/optim/ug/linprog.html)
 *
 *         min   -1 * x - 1/3 *  y
 *         s.t.       x +        y <=  2
 *                    x + 0.25 * y <=  1
 *                    x -        y <=  2
 *            -0.25 * x -        y <=  1
 *                   -x -        y <= -1
 *                   -x +        y <=  2
 *
 * The solution is (2/3, 4/3) and the value is -10/9
 */
#define m 6        /* number of constraints */
#define n 2        /* number of variables   */

double obj[] = {-1., -1./3., 0.5};
double constraint_1_coef[] = {  1.,  1.   };
double constraint_2_coef[] = {  1.,  0.25 };
double constraint_3_coef[] = {  1., -1.   };
double constraint_4_coef[] = { -0.25, -1.   };
double constraint_5_coef[] = { -1., -1.   };
double constraint_6_coef[] = { -1.,  1.   };

struct impf_LinearConstraint constraints[] = {
	{ "", constraint_1_coef,  2., impf_CONS_T_LE },
	{ "", constraint_2_coef,  1., impf_CONS_T_LE },
	{ "", constraint_3_coef,  2., impf_CONS_T_LE },
	{ "", constraint_4_coef,  1., impf_CONS_T_LE },
	{ "", constraint_5_coef, -1., impf_CONS_T_LE },
	{ "", constraint_6_coef,  2., impf_CONS_T_LE }
};
struct impf_VariableBound bounds[] = {
	{ "x", __impf_NINF__, __impf_INF__, impf_BOUND_T_FR, impf_VAR_T_REAL },
	{ "y", __impf_NINF__, __impf_INF__, impf_BOUND_T_FR, impf_VAR_T_REAL }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = impf_lp_simplex(obj, constraints, bounds, m, n, "bland", 1000, x, &value, &code);

	printf("error = %i\n", code);
	assert(state == EXIT_SUCCESS);
	assert(__impf_ABS__(value + 10. / 9.) < 1e-8);
	assert(__impf_ABS__(x[0] - 2./3.) < 1e-8);
	assert(__impf_ABS__(x[1] - 4./3.) < 1e-8);
	return 0;
}
