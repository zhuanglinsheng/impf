#include <stdio.h>
#include <stdlib.h>
#include <impf/impf.h>

int main(void)
{
    /* matrix A of 3 by 3 */
    double A[3 * 3] = {1.5, 3.0, 2.2, 4.5, 2.1, 6.0, 0.9, 2.1, 3.3};

    /* vector b of 3 */
    double b[3] = {1.0, 1.0, 1.0};

    /* solve the linear equations of 3 unknowns */
    double x[3];
    if(IMPF_EXIT_SUCCESS != impf_subrt_lsolve_3A3_b3(
        A[0], A[1], A[2],
        A[3], A[4], A[5],
        A[6], A[7], A[8],
        b[0], b[1], b[2], x))
    {
        exit(-1);
    }

    /* calculate the residule */
    double err[3] = {b[0], b[1], b[2]}; /* err = b */
    impf_t_matrix matA = {A, 3, 3, IMPF_MAT_ROW_MAJOR};
    impf_t_vector vecx = {x, 3};
    impf_t_vector vecerr = {err, 3};
    impf_subrt_mv_mul(1, &matA, NULL, 'N', &vecx, 1, -1, &vecerr, 1); /* Ax - b */

    /* display the result. */
    FILE * out = fopen("eg_linalg_lsolve_3A3_3b.out", "w");
    fprintf(out, "The solution of 'Ax = b' is (%f, %f, %f)\n", x[0], x[1], x[2]);
    fprintf(out, "Ax - b = \n");
    impf_fprint_vec_sci(out, &vecerr);
    fclose(out);
    return 0;
}
