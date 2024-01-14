/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#ifndef __IMPF_BASIC_H__
#define __IMPF_BASIC_H__

#ifdef __cpluscplus
extern "C" {
#endif /* __cplusplus */

#include <assert.h>
#include <stddef.h>

#define __impf_ABS__(x) ((x) >= 0 ? (x) : (-(x)))
#define __impf_MAX__(x, y) ((x) >= (y) ? (x) : (y))
#define __impf_MIN__(x, y) ((x) <= (y) ? (x) : (y))

#define __impf_INF__ (1. / 0.)
#define __impf_NINF__ (-1. / 0.)

#define impf_Success			0
#define impf_MemoryAllocError		1
#define impf_CondUnsatisfied		2
#define impf_ExceedIterLimit		3
#define impf_Singularity		4
#define impf_OverDetermination		5
#define impf_Unboundedness		6
#define impf_Infeasibility		7
#define impf_Degeneracy			8
#define impf_PrecisionError		9

/*******************************************************************************
 * Utils of "arri-info-family"
 ******************************************************************************/

int is_in_arri(const int idx, const int *idxset, const int len);
int maxabs_arri(const int *arr, const int len, const int inc);

/*******************************************************************************
 * Utils of "arrd-info-family"
 *
 * Note: there is no float error concern for
 *     1. taking absolute value (just change sign)
 *     2. finding the max value or index (just compare values)
 ******************************************************************************/

int argmaxabs_arrd(const double *arr, const int len, const int inc);
double maxabs_arrd(const double *arr, const int len, const int inc);
double maxabs_arrd_gap(const double *arr1, const double *arr2, const int len, const int inc);

/*******************************************************************************
 * Utils of "prt-family"
 ******************************************************************************/

void impf_prt_arri(const int *arr, const int len, const int inc);
void impf_prt_arrl(const long *arr, const int len, const int inc);
void impf_prt_arrd(const double *arr, const int len, const int inc, const int sci);
void impf_prt_arrld(const long double *arr, const int len, const int inc, const int sci);

void impf_prt_matd(const double *mat, const int ld, const int nrow, const int ncol);

/*******************************************************************************
 * Utils of "memory-family"
 ******************************************************************************/

void *impf_malloc(size_t size);
void impf_free(void *ptr);
void *impf_memset(void *str, int c, size_t n);
void *impf_memcpy(void *dest, const void *src, size_t n);
int impf_memcmp(const void *str1, const void *str2, size_t n);
size_t impf_strcspn(const char *str1, const char *str2);
size_t impf_strlen(const char *str);
double impf_atof(const char* str);

/*******************************************************************************
 * Utils of "state-family"
 ******************************************************************************/

#define impf_EXIT_SUCCESS 0
#define impf_EXIT_FAILURE 1

/*******************************************************************************
 * Utils of "mode-family"
 ******************************************************************************/

/* developing mode (only for development) */
/* #define IMPF_MODE_DEV */

/* debug mode (only for for development) */
/* #define IMPF_MODE_DEBUG */

#ifdef __cpluscplus
}
#endif /* __cpluscplus */

#endif
