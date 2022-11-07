#include <stdio.h>
#include <impf/impf.h>

/* m = 4, n = 4 */
double m_4_4_data[4 * 4] = {
    1., 2., 3., 4.,
    5., 4., 3., 2.,
    9., 8., 7., 6.,
    5., 6., 7., 8.,
};
double m_4_4_data_t[4 * 4] = {
    1., 5., 9., 5.,
    2., 4., 8., 6.,
    3., 3., 7., 7.,
    4., 1., 6., 8.,
};
double v_4_data[4] = {
    10., 14., 30., 26.,
};
impf_t_matrix mat_4_4 = {m_4_4_data_t, 4, 4, IMPF_MAT_COL_MAJOR};
impf_t_vector vec_4 = {v_4_data, 4};


int main(void)
{
    return 0;
}
