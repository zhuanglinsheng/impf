/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/lp.h>
#include <stdio.h>

#define m 27        /* number of constraints */
#define n 32        /* number of variables   */

double obj[] = { 0., -0.4, 0., 0., 0., 0., 0., 0.,
		 0.,  0., 0., 0., -0.32, 0., 0., 0.,
		-0.6,  0., 0., 0., 0., 0., 0., 0.,
		 0.,  0., 0., 0., -0.48, 0., 0., 10.};
double cons01_coef[] = {-1., 1., 1., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons02_coef[] = {-1.06, 0., 0., 1., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons03_coef[] = { 1., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons04_coef[] = { 0., -1., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 1.40, 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons05_coef[] = { 0., 0., 0., 0., -1., -1., -1., -1.,
			 0., 0., 0., 0., 1., 1., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons06_coef[] = { 0., 0., 0., 0., -1.06, -1.06, -0.96, -0.86,
			 0., 0., 0., 0., 0., 0., 1., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons07_coef[] = { 0., 0., 0., 0., 1., 0., 0., 0.,
			-1., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons08_coef[] = { 0., 0., 0., 0., 0., 1., 0., 0.,
			 0., -1., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons09_coef[] = { 0., 0., 0., 0., 0., 0., 1., 0.,
			 0., 0., -1., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons10_coef[] = { 0., 0., 0., 0., 0., 0., 0., 1.,
			 0., 0., 0., -1., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons11_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., -1.,
			 1., 1., 1., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons12_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., -0.43,
			 0., 0., 0., 1., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons13_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 1.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons14_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			-1., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 1.40, 0., 0., 0.};
double cons15_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., -0.43, -0.43, -0.39, -0.37,
			 0., 0., 0., 0., 0., 0., 1., 0.};
double cons16_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 1., 1., 1., 1.,
			 0., 0., 0., 0., -1., 1., 0., 1.};
double cons17_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 1., 0., 0., 0.,
			-1., 0., 0., 0., 0., 0., 0., 0.};
double cons18_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 1., 0., 0.,
			 0., -1., 0., 0., 0., 0., 0., 0.};
double cons19_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 1., 0.,
			 0., 0., -1., 0., 0., 0., 0., 0.};
double cons20_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 1.,
			 0., 0., 0., -1., 0., 0., 0., 0.};
double cons21_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 2.364, 2.386, 2.408, 2.429, 0., 0., 0., 0.,
			 0., 0., -1., 0., 0., 0., 0., 0.,
			 2.191, 2.219, 2.249, 2.279, 0., 0., 0., 0.};
double cons22_coef[] = { 0., 0., -1., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.109,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons23_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., -1., 0., 0.,
			 0., 0., 0., 0., 0.109, 0.108, 0.108, 0.107,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons24_coef[] = { 0.301, 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., -1., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons25_coef[] = { 0., 0., 0., 0., 0.301, 0.313, 0.313, 0.326,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., -1., 0., 0.};
double cons26_coef[] = { 0., 0., 0., 1., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 1., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.};
double cons27_coef[] = { 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 1., 0.,
			 0., 0., 0., 0., 0., 0., 0., 0.,
			 0., 0., 0., 0., 0., 0., 1., 0.};
struct impf_LinearConstraint constraints[] = {
	{ "", cons01_coef,   0., impf_EQ }, { "", cons02_coef,   0., impf_EQ },
	{ "", cons03_coef,  80., impf_LE }, { "", cons04_coef,   0., impf_LE },
	{ "", cons05_coef,   0., impf_EQ }, { "", cons06_coef,   0., impf_EQ },
	{ "", cons07_coef,  80., impf_LE }, { "", cons08_coef,   0., impf_LE },
	{ "", cons09_coef,   0., impf_LE }, { "", cons10_coef,   0., impf_LE },
	{ "", cons11_coef,   0., impf_EQ }, { "", cons12_coef,   0., impf_EQ },
	{ "", cons13_coef, 500., impf_LE }, { "", cons14_coef,   0., impf_LE },
	{ "", cons15_coef,   0., impf_EQ }, { "", cons16_coef,   0., impf_EQ },
	{ "", cons17_coef,  44., impf_LE }, { "", cons18_coef, 500., impf_LE },
	{ "", cons19_coef,   0., impf_LE }, { "", cons20_coef,   0., impf_LE },
	{ "", cons21_coef,   0., impf_LE }, { "", cons22_coef,   0., impf_LE },
	{ "", cons23_coef,   0., impf_LE }, { "", cons24_coef,   0., impf_LE },
	{ "", cons25_coef,   0., impf_LE }, { "", cons26_coef, 310., impf_LE },
	{ "", cons27_coef, 300., impf_LE }
};

int main(void)
{
	/* call simplex subroutine */
	double x[n], value;
	int code;
	int state = impf_lp_simplex(obj, constraints, NULL, m, n, "bland", 1000, x, &value, &code);

	printf("error code = %u\n", code);
	assert(state == EXIT_SUCCESS);
	printf("value = %.10e\nSolution = ", value);
	impf_prt_arrd(x, n, 1, 0);
	printf("\n");
	return 0;
}
