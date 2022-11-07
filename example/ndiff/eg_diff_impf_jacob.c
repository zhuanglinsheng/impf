#include <stdio.h>
#include <math.h>
#include <impf/impf.h>

void mfn(const unsigned int m, const unsigned int n, const double * x, double * y)
{
    unsigned int i, j;

    for(i = 0; i<m; i++)
    {
        y[i] = 0;

        for(j = 0; j<n; j++)
        {
            y[i] += pow(x[j] - i - j, 2);
        }
    }
}

void jac_mfn(const unsigned int m, const unsigned int n, const double * x, double * jac)
{
    unsigned int i, j;

    for(i = 0; i<m; i++)
    {
        for(j = 0; j<n; j++)
        {
            jac[i * n + j] = 2 * (x[j] - i - j);
        }
    }
}

int main(void)
{
    FILE * out = fopen("eg_diff_impf_jacob.out", "w");
    if(NULL == out) return -1;

    /* points: m = 4, n = 5 */
    double point[5] = {1./1, 1.0/2, 1.0/3, 1.0/4, 1.0/5};
    double buf[5 + 4 * 4];
    double njac_data[4 * 5]; /* numerical result  */
    impf_t_matrix njac = {njac_data, 4, 5, IMPF_MAT_ROW_MAJOR};
    double ajac_data[4 * 5]; /* analytical result */
    impf_t_matrix ajac = {ajac_data, 4, 5, IMPF_MAT_ROW_MAJOR};

    /* Analytical version */
    jac_mfn(4, 5, point, ajac_data);
    fprintf(out, "Analytical:\n");
    impf_fprint_mat_sci(out, &ajac, 100, 0);

    /* Numerical version */
    impf_subrt_jacob_mfn(mfn, 4, point, 5, buf, &njac);
    fprintf(out, "Numerical:\n");
    impf_fprint_mat_sci(out, &njac, 100, 0);

    fclose(out);
    return 0;
}
