/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/fmin_lp.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

void netlib_bm(const char *bm, const int n, const double fval)
{
	clock_t start, stop;
	double usec = 0, msec = 0, sec = 0;

	double *x, value;
	int code;
	int state;
	struct impf_Model_LP *model = NULL;
	char file[80];

	memset(file, '\0', 80);
	memcpy(file, "../../include/netlib/lp/data/", 29);
	memcpy(file + 29, bm, n);
	model = impf_lp_readmps(file);
	assert(model != NULL);
	printf("\nBenchmark = \"%s\"\n", bm);
	printf("m = %i, n = %i\n", model->m, model->n);
	x = malloc((model->n) * sizeof(double));

	start = clock();
	state = impf_lp_simplex_wrp(model, "bland", 2000, x, &value, &code);
	stop = clock();
	usec = (double)((stop - start) / (CLOCKS_PER_SEC / 1000000));
	msec = (double)((stop - start) / (CLOCKS_PER_SEC / 1000));
	sec = (double)((stop - start) / CLOCKS_PER_SEC);

	printf("Duration = %.fus, %.4fms, %.4fs\n", usec, msec, sec);
	printf("Error code = %u\n", code);
	assert(state == EXIT_SUCCESS);
	printf("value = %.11e\n\n", value);
	assert(value <= fval);
	impf_lp_free(model);
	free(x);
}

int main(void)
{
	/* degenerate
	netlib_bm("25fv47", 6, 5.5018458883E+03);
	netlib_bm("80bau3b", 7, 9.8723216072E+05);
	*/

	/* Low precision
	netlib_bm("adlittle", 8, 2.2549496316E+05);
	*/

	/* ok
	netlib_bm("afiro", 5, -4.6475314286E+02);
	netlib_bm("agg", 3, -3.5991767287E+07);
	netlib_bm("agg2", 4, -2.0239252356E+07);
	netlib_bm("agg3", 4, 1.0312115935E+07);
	*/

	/* infeasible or unbounded
	netlib_bm("bandm", 5, -1.5862801845E+02);
	netlib_bm("beaconfd", 8, 3.3592485807E+04);
	netlib_bm("blend", 5, -3.0812149846E+01);
	*/

	/* degenerate
	netlib_bm("bnl1", 4, 1.9776292856E+03);
	netlib_bm("bnl2", 4, 1.8112365404E+03);
	netlib_bm("brandy", 6, 1.5185098965E+03);
	*/

	/* with BOUNDS, RANGES, empty RHS
	netlib_bm("boeing1", 7, -3.3521356751E+02);
	netlib_bm("boeing2", 7, -3.1501872802E+02);
	netlib_bm("bore3d", 6, 1.3730803942E+03);
	*/
	return 0;
}
