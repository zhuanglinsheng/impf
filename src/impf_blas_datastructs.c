/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <impf/impf_basic.h>


void
impf_mat_zeros(impf_t_matrix * mat)
{
    memset(mat->data, 0, sizeof(double) * mat->nrow * mat->ncol);
}

void
impf_mat_ones(impf_t_matrix * mat)
{
    memset(mat->data, 0, sizeof(double) * mat->nrow * mat->ncol);
    unsigned long int i;

    for(i = 0; i<mat->ncol * mat->nrow; i++)
    {
        *(mat->data + i) = 1;
    }
}

void
impf_mat_eyes(impf_t_matrix * mat)
{
    unsigned int i, LD, rank;
    
    i = 0;
    LD = IMPF_MAT_COL_MAJOR == mat->major ? mat->nrow : mat->ncol;
    rank = mat->nrow <= mat->ncol ? mat->nrow : mat->ncol;
    
ASSIGNMENT:
    *(mat->data + i * LD + i) = 1;
    i ++;
    if(i < rank) goto ASSIGNMENT;
}

void
impf_mat_transpose(impf_t_matrix * mat)
{
    unsigned int tmp = mat->ncol;
    mat->ncol = mat->nrow;
    mat->nrow = tmp;
    mat->major = IMPF_MAT_COL_MAJOR == mat->major ? \
        IMPF_MAT_ROW_MAJOR : IMPF_MAT_COL_MAJOR;
}

void
impf_mat_transmajor(impf_t_matrix * mat)
{
    unsigned long int size = sizeof(double) * mat->ncol * mat->nrow;
    double * data = malloc(size);
    int i, j;

    if(NULL == data)
    {
        printf("Memory allocation failure in the function 'impf_mat_transmajor(impf_t_matrix *)'.\n");
        exit(EXIT_FAILURE);
    }
    for(i = 0; i<mat->nrow; i++)
    {
        for(j = 0; j<mat->ncol; j++)
        {
            *(data + j * mat->nrow + i) = *(mat->data + i * mat->ncol + j);
        }
    }
    memcpy(mat->data, data, size);
    free(data);
    mat->major = IMPF_MAT_COL_MAJOR == mat->major ? \
        IMPF_MAT_ROW_MAJOR : IMPF_MAT_COL_MAJOR;
}

/* just for testing */
void
impf_mat_transmajor_inplace(impf_t_matrix * mat)
{
    if(1 >= mat->nrow || 1 >= mat->ncol)
    {
        mat->major = 'T' == mat->major ? IMPF_MAT_ROW_MAJOR : IMPF_MAT_COL_MAJOR;
        return;
    }

    if(IMPF_MAT_ROW_MAJOR == mat->major)
    {
        if(mat->ncol == mat->nrow) /* for squared matrix */
        {
            double tmp;
            unsigned int i, j;

            for(i = 0; i < mat->nrow; i++)
            {
                for(j = 0; j < i; j++)
                {
                    tmp = *(mat->data + i * mat->ncol + j);
                    *(mat->data + i * mat->ncol + j) = *(mat->data + j * mat->nrow + i);
                    *(mat->data + j * mat->nrow + i) = tmp;
                }
            }
        }
        else /* for non-squared matrix */
        {
            double vstart;
            unsigned long int a, pa, start, lenm1;
            unsigned long int len = mat->nrow * mat->ncol;
            char * recorder = malloc(len);

            if(NULL == recorder)
            {
                printf("Memory allocation failure in the function 'impf_mat_transmajor_2(impf_t_matrix *)'.\n");
                exit(EXIT_FAILURE);
            }

            memset(recorder, 0, len);
            start = 1;
            lenm1 = len - 1;
            recorder[0] = 1;
            recorder[lenm1] = 1;
        LOOP:
            for(; start<lenm1 && recorder[start]; start++);
            if(start == lenm1) goto END;
            vstart = *(mat->data + start);
            a = start;
            pa = (mat->ncol * a) % lenm1;

            while(start != pa)
            {
                *(mat->data + a) = *(mat->data + pa);
                recorder[a] = 1;
                a = pa;
                pa = (mat->ncol * a) % lenm1;
            }
            *(mat->data + a) = vstart;
            recorder[a] = 1;
            goto LOOP;
        END:
            free(recorder);
        }
        mat->major = IMPF_MAT_COL_MAJOR;
    }
    else
    {
        impf_mat_transpose(mat);
        impf_mat_transmajor(mat);
        impf_mat_transpose(mat);
    }
}

void
impf_submat_subtract(
    const impf_t_matrix * mat, const impf_t_submatinfo * submat, const char trans,
    int * M, int * N, double ** data)
{
    assert(NULL != mat);
    assert(NULL != M);
    assert(NULL != N);
    assert(NULL != data);

    int row_room, col_room;

    if(NULL != submat)
    {
        row_room = mat->nrow - submat->iloc;
        row_room = row_room <= submat->nrow ? row_room : submat->nrow;
        col_room = mat->ncol - submat->jloc;
        col_room = col_room <= submat->ncol ? col_room : submat->ncol;
        
        if(IMPF_MAT_COL_MAJOR == mat->major)
        {
            *data = mat->data + submat->jloc * mat->nrow + submat->iloc;
        }
        else
        {
            *data = mat->data + submat->iloc * mat->ncol + submat->jloc;
        }
    }
    else
    {
        row_room = mat->nrow;
        col_room = mat->ncol;
        *data = mat->data;
    }

    if(0 >= row_room || 0 >= col_room)
    {
        *M = 0;
        *N = 0;
        *data = NULL;
        return;
    }

    if('N' == trans)
    {
        *M = mat->nrow <= row_room ? mat->nrow : row_room;
        *N = mat->ncol <= col_room ? mat->ncol : col_room;
    }
    else
    {
        *N = mat->nrow <= row_room ? mat->nrow : row_room;
        *M = mat->ncol <= col_room ? mat->ncol : col_room;
    }
}

void 
impf_print_vec_sci(const impf_t_vector * vec)
{
    impf_fprint_vec_sci(stdout, vec);
}

void
impf_fprint_vec_sci(FILE * out, const impf_t_vector * vec)
{
    assert(NULL != out);
    assert(NULL != vec);

    unsigned int i;
    double * data = vec->data;

    for(i = 0; i<vec->dim; i++)
    {
        fprintf(out, "%16.9e\n", *data);
        data += 1;
    }
}


void 
impf_print_mat_sci(const impf_t_matrix * mat, const unsigned int width, const unsigned foot)
{
    impf_fprint_mat_sci(stdout, mat, width, foot);
}

void
impf_fprint_mat_sci(FILE * out, const impf_t_matrix * mat, const unsigned int width, const unsigned foot)
{
    assert(NULL != out);
    assert(NULL != mat);
    assert(width >= 40);

    unsigned int print_num = (width - 4) / 18;
    unsigned int i, j, base;
    double * data = mat->data;

    if(IMPF_MAT_ROW_MAJOR == mat->major)
    {
        for(i = 0; i<mat->nrow; i++)
        {
            base = i * mat->ncol;

            for(j = 0; j<mat->ncol; j++)
            {
                if(j == print_num - 1 && j < mat->ncol - 1)
                {
                    fprintf(out, "..  %16.9e", data[base + mat->ncol - 1]);
                    break;
                }
                else
                {
                    fprintf(out, "%16.9e  ", data[base + j]);
                }
            }
            fprintf(out, "\n");
        }
    }
    else
    {
        base = (mat->ncol - 1) * (mat->nrow);

        for(i = 0; i<mat->nrow; i++) /* Rows */
        {
            for(j = 0; j<mat->ncol; j++) /* Cols */
            {
                if(j == print_num - 1 && j < mat->ncol - 1)
                {
                    fprintf(out, ".. %16.9e", data[base + i]);
                    break;
                }
                else
                {
                    fprintf(out, "%16.9e  ", data[j * mat->nrow + i]);
                }
            }
            fprintf(out, "\n");
        }
    }
    if(foot)
    {
        fprintf(out, "Matrix of (%ux%u)", mat->nrow, mat->ncol);
        if(IMPF_MAT_ROW_MAJOR == mat->major)
            fprintf(out, ", row-major\n");
        else
            fprintf(out, ", col-major\n");
    }
}
