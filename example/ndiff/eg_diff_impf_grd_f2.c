#include <stdio.h>
#include <impf/impf.h>

/* 2-dimensional real function f(x,y) = x^2 + 2xy + y^2.
 * Match the real function type `impf_t_f2`. 
 */
double f2(const double x, const double y)
{
    return x * x + 2 * x * y + y * y;
}

int main(void)
{
    FILE * out = fopen("eg_diff_impf_grd_f2.out", "w");
    if(NULL == out) return -1;

    double x = 1.0 / 3.0;
    double y = 5.0 / 7.0;
    double grd [2];
    impf_subrt_gradient_f2(f2, x, y, grd);
    fprintf(out, "The gradient of f2 at (%f,%f) is (%f,%f)\n", x, y, grd[0], grd[1]);

    fclose(out);
    return 0;
}
