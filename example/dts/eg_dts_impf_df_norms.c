#include <stdio.h>
#include <impf/impf.h>

double x[10] = {1.0, 3.2, 5.5, -0.9, -7.8, 2.3, 7.0, 3.3, 4.0, 5.1};
impf_t_vector vec = {x, 10};

int main(void)
{
    FILE * out = fopen("eg_dts_impf_df_norms.out", "w");
    if(NULL == out) return -1;

    fprintf(out, "|v|_1   = %f\n", impf_df_norm1(&vec, 1));
    fprintf(out, "|v|_2   = %f\n", impf_df_norm2(&vec, 1));
    fprintf(out, "|v|_inf = %f\n", impf_df_norminf(&vec, 1));

    fclose(out);
    return 0;
}
