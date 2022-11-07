#include <stdio.h>
#include <impf/impf.h>

double x[10] = {1.0, 3.2, 5.5, -0.9, -7.8, 2.3, 7.0, 3.3, 4.0, 5.1};
double y[10] = {4.4, 7.1, 1.2, 3.7, 8.0, -2.0, 6.1, -3.5, 9.1, 4.7};
impf_t_vector vx = {x, 10};  /* x is a 10-dimensional double array */
impf_t_vector vy = {y, 10};  /* y is a 10-dimensional double array */

int main(void)
{
    double buf[10];
    FILE * out = fopen("eg_dts_impf_df_distances.out", "w");
    if(NULL == out) return -1;

    fprintf(out, "|x[0   ] - y[0   ]| = %f\n", impf_df_abs_sub(*x, *y));
    fprintf(out, "|x[0..1] - y[0..1]| = %f\n", impf_df_distance_norminf_2d(x, y));
    fprintf(out, "|x[0..2] - y[0..2]| = %f\n", impf_df_distance_norminf_3d(x, y));
    fprintf(out, "|v1 - v2|_1         = %f\n", impf_df_distance_norm1(&vx, &vy, 1, 1, 10, buf));
    fprintf(out, "|v1 - v2|_2         = %f\n", impf_df_distance_norm2(&vx, &vy, 1, 1, 10, buf));
    fprintf(out, "|v1 - v2|_inf       = %f\n", impf_df_distance_norminf(&vx, &vy, 1, 1, 10, buf));

    fclose(out);
    return 0;
}
