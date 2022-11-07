#include <impf/impf.h>

/* Functional F(y(x), x) = y + x - y * x + x^2 + y^2
 * Matching the wrapper type `impf_t_1f1_wrapper`.
 */
void wrp_1f1(const impf_t_f1 f,  /* the implicit function f */
             const double x,     /* the variant x */
             double *res)        /* the evaluation result of F */
{
    double y = f(x);
    *res = y + x - y * x + x*x + y * y;
}
