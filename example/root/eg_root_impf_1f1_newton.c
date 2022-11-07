#include <stdio.h>
#include <math.h>
#include <impf/impf.h>

double f(const double x)
{
    return 1 + x + pow(1 + x, 2) + pow(1 + x, 3);
}

int main(void)
{
    FILE * out = fopen("eg_root_impf_1f1_newton.out", "w");
    if(NULL == out) return -1;

    double x = -0.0; /* the initial guess */

    if(IMPF_EXIT_SUCCESS == impf_subrt_1f1_root_newton(f, &x, 1e-8, 1e-8, 1000))
    {
        fprintf(out, "x    = %f\n", x);
        fprintf(out, "f(x) = %f\n", f(x));
    }
    else
    {
        fprintf(out, "Newton method fail to converge.\n");
    }
    fclose(out);
    return 0;
}
