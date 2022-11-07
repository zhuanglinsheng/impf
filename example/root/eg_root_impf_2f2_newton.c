#include <stdio.h>
#include <math.h>
#include <impf/impf.h>

double f1(const double x, const double y)
{
    return 1.5 + x + y + pow(0.5 + x, 2) + pow(0.5 + y, 3);
}
double f2(const double x, const double y)
{
    return 2.0 + x + y + pow(1.0 + y, 2) + pow(1.0 + x, 3);
}

int main(void)
{
    FILE * out = fopen("eg_root_impf_2f2_newton.out", "w");
    if(NULL == out) return -1;

    double x = -0.0; /* the initial guess */
    double y = -0.0;

    if(IMPF_EXIT_SUCCESS == impf_subrt_2f2_root_newton(f1, f2, &x, &y, 1e-8, 1e-8, 1000))
    {
        fprintf(out, "x = %f, y = %f\n", x, y);
        fprintf(out, "f1(x, y) = %f\n", f1(x, y));
        fprintf(out, "f2(x, y) = %f\n", f2(x, y));
    }
    else
    {
        fprintf(out, "Newton method fail to converge.\n");
    }
    fclose(out);
    return 0;
}
