#include <stdio.h>
#include <stdlib.h>
#include <impf/impf.h>

int main(void)
{
    /* matrix A of 2 by 2 */
    double A[2 * 2] = {1.5, 3.0, 4.5, 2.2};

    /* vector b of 2 */
    double b[2] = {1.0, 1.0};

    /* solve the linear equations of 2 unknowns */
    double x[2];
    if(IMPF_EXIT_SUCCESS != impf_subrt_lsolve_2A2_b2(A[0], A[1], A[2], A[3], b[0], b[1], x))
    {
        exit(-1);
    }

    /* calculate the residule */
    double err[2] = {b[0], b[1]}; /* err = b */
    impf_t_matrix matA = {A, 2, 2, IMPF_MAT_ROW_MAJOR};
    impf_t_vector vecx = {x, 2};
    impf_t_vector vecerr = {err, 2};
    impf_subrt_mv_mul(1, &matA, NULL, 'N', &vecx, 1, -1, &vecerr, 1); /* Ax - b */

    /* display the result. */
    FILE * out = fopen("eg_linalg_lsolve_2A2_2b.out", "w");
    fprintf(out, "The solution of 'Ax = b' is (%f, %f)\n", x[0], x[1]);
    fprintf(out, "Ax - b = \n");
    impf_fprint_vec_sci(out, &vecerr);
    fclose(out);
    return 0;
}
