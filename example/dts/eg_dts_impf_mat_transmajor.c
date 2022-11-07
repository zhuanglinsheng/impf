#include <stdio.h>
#include <impf/impf.h>

int
main(void)
{
    /* Create a matrix */
    double data[] = {1., 2., 3., 4., 5., 6.};
    impf_t_matrix mat = {data, 2, 3, IMPF_MAT_ROW_MAJOR};

    /* Change the major of a matrix */
    impf_mat_transmajor(&mat);

    /* Print the matrix */
    FILE * out = fopen("eg_dts_impf_mat_transmajor.out", "w");
    if(NULL == out) return -1;
    impf_fprint_mat_sci(out, &mat, 100, 1);
    fclose(out);
    return 0;
}
