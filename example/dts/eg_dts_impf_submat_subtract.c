#include <stdio.h>
#include <impf/impf.h>

int
main(void)
{
    FILE * out = fopen("eg_dts_impf_submat_subtract.out", "w");
    if(NULL == out) return -1;

    /* Create a matrix */
    double data[] = {1., 2., 3., 4., 5., 6., 7., 8., 9.};
    impf_t_matrix mat = {data, 3, 3, IMPF_MAT_ROW_MAJOR};

    /* Print the matrix */
    fprintf(out, "Outer matrix:\n");
    impf_fprint_mat_sci(out, &mat, 100, 1);

    /* Declare a submatrix: starting from (0,1) of dimension 3 by 2. */
    impf_t_submatinfo submat = {0, 1, 3, 2};

    /* Subtract a submatrix from matrix */
    double * subdata;
    int M, N;
    impf_submat_subtract(&mat, &submat, 'T', &M, &N, &subdata);

    /* Print the matrix */
    fprintf(out, "\nSubmatrix information:\n");
    fprintf(out, "M = %i, N = %i, First Element = %f\n", M, N, *subdata);

    fclose(out);
    return 0;
}
