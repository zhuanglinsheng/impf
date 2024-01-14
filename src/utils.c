/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/utils.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int is_in_arri(const int idx, const int *idxset, const int len)
{
	int i;

	for (i = 0; i < len; i++) {
		if (idx == idxset[i])
			return 1;
	}
	return 0;
}

int maxabs_arri(const int *arr, const int len, const int inc)
{
	int i, ele, maxv = __impf_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i++) {
		ele = __impf_ABS__(arr[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

int argmaxabs_arrd(const double *arr, const int len, const int inc)
{
	int i, idx = 0;
	double ele, maxv = __impf_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __impf_ABS__(arr[i]);
		if (ele > maxv) {
			maxv = ele;
			idx = i;
		}
	}
	return idx;
}

double maxabs_arrd(const double *arr, const int len, const int inc)
{
	int i;
	double ele, maxv = __impf_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __impf_ABS__(arr[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

double maxabs_arrd_gap(const double *arr1, const double *arr2, const int len, const int inc)
{
	int i;
	double ele, maxv = __impf_NINF__;

	assert(inc >= 1);

	for (i = 0; i < len; i += inc) {
		ele = __impf_ABS__(arr1[i] - arr2[i]);
		if (ele > maxv)
			maxv = ele;
	}
	return maxv;
}

void impf_prt_arri(const int *arr, const int len, const int inc)
{
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		printf("%i", arr[j]);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void impf_prt_arrl(const long *arr, const int len, const int inc)
{
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		printf("%li", arr[j]);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void impf_prt_arrd(const double *arr, const int len, const int inc, const int sci)
{
	int j;
	double ele;

	printf("[");
	for (j = 0; j < len; j += inc) {
		ele = arr[j];
		if (-0. == ele)
			ele = 0.;
		if (sci)
			printf("%e", ele);
		else
			printf("%f", ele);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void impf_prt_arrld(const long double *arr, const int len, const int inc, const int sci)
{
	long double ele;
	int j;

	printf("[");
	for (j = 0; j < len; j += inc) {
		ele = arr[j];
		if (-0. == ele)
			ele = 0.;
		if (sci)
			printf("%Le", ele);
		else
			printf("%Lf", ele);
		if (j + inc < len)
			printf(", ");
	}
	printf("]");
}

void impf_prt_matd(const double *mat, const int ld, const int nrow, const int ncol)
{
	int i, j;

	for (i = 0; i < nrow; i++) {
		for (j = 0; j < ncol; j++) {
			printf("%e", mat[j + i * ld]);
			if (j < ncol - 1)
				printf(", ");
		}
		if(i < nrow - 1)
			printf("\n");
	}
}

void *impf_malloc(size_t size)
{
	return malloc(size);
}

void impf_free(void *ptr)
{
	free(ptr);
}

void *impf_memset(void *str, int c, size_t n)
{
	return memset(str, c, n);
}

void *impf_memcpy(void *dest, const void *src, size_t n)
{
	return memcpy(dest, src, n);
}

int impf_memcmp(const void *str1, const void *str2, size_t n)
{
	return memcmp(str1, str2, n);
}

size_t impf_strcspn(const char *str1, const char *str2)
{
	return strcspn(str1, str2);
}

size_t impf_strlen(const char *str)
{
	return strlen(str);
}

double impf_atof(const char* str)
{
	return atof(str);
}
