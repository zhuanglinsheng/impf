#include <stdio.h>
#include <impf/impf.h>

int main(void)
{
    /* Create a matrix of dimensions 2 by 3. */
    double data[] = {1., 2., 3., 4., 5., 6.};
    impf_t_matrix mat = {data, 2, 3, IMPF_MAT_ROW_MAJOR};

    /* Display the matrix. */
    FILE * out = fopen("eg_dts_impf_t_mat.out", "w");
    if(NULL == out) return -1;
    impf_fprint_mat_sci(out, &mat, 100, 1);
    fclose(out);
    return 0;
}
