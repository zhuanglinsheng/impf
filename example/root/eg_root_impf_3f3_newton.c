#include <stdio.h>
#include <math.h>
#include <impf/impf.h>

double f1(const double x, const double y, const double z)
{
    return 1.5 + x + y + pow(0.5 + x, 2) + pow(0.5 + y, 3);
}
double f2(const double x, const double y, const double z)
{
    return 2.0 + x + y + pow(1.0 + y, 2) + pow(1.0 + x, 3);
}
double f3(const double x, const double y, const double z)
{
    return 2 + x + y + z;
}

int main(void)
{
    FILE * out = fopen("eg_root_impf_3f3_newton.out", "w");
    if(NULL == out) return -1;

    double x = -0.0; /* the initial guess */
    double y = -0.0;
    double z = -0.0;

    if(IMPF_EXIT_SUCCESS == impf_subrt_3f3_root_newton(f1, f2, f3, &x, &y, &z, 1e-8, 1e-8, 1000))
    {
        fprintf(out, "x = %f, y = %f, z = %f\n", x, y, z);
        fprintf(out, "f1(x, y, z) = %f\n", f1(x, y, z));
        fprintf(out, "f2(x, y, z) = %f\n", f2(x, y, z));
        fprintf(out, "f3(x, y, z) = %f\n", f3(x, y, z));
    }
    else
    {
        fprintf(out, "Newton method fail to converge.\n");
    }
    fclose(out);
    return 0;
}
