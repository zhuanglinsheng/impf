#include <impf/impf.h>

/* The definition of 2-dimensional heat equation, following the wrapper  
 * type `impf_t_1f2_wrapper`.
 * 
 * Suppose \f$u(t,x):\mathbb{R}^2\to\mathbb{R}\f$, the classical heat 
 *     equation of 2-dimensional is 
 * \f[
 * \frac{u(t+h, x) - u(t-h, x)}{2h} = \frac{u(t,x-h)-2u(t,x)+u(t,x+h)}{h^2}
 * \f]
 */
void heat2D(const impf_t_f2 u,  /* the implicit function u */
            const double t,     /* time */
            const double x,     /* location on the line */
            double *res)        /* the residule of the equation */
{
    double h = 1e-6;
    *res = (u(t + h, x) - u(t - h, x)) / (2 * h) - \
           (u(t, x - h) - 2 * u(t, x) + u(t, x + h)) / (h * h);
}

int main(void)
{
    return 0;
}
