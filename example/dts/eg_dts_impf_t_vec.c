#include <stdio.h>
#include <impf/impf.h>

int main(void)
{
    /* Create a vector of dimension 6 */
    double data[] = {1., 2., 3., 4., 5., 6.};
    impf_t_vector vec = {data, 6};

    /* Display the vector */
    FILE * out = fopen("eg_dts_impf_t_vec.out", "w");
    if(NULL == out) return -1;
    impf_fprint_vec_sci(out, &vec);
    fclose(out);
    return 0;
}
