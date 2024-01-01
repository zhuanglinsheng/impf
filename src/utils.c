/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/utils.h>
#include <stdio.h>

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

/* Find the maximum absolute value of a row
 */
long double table_abs_rowmax(const long double *table, const int ldtable,
			     const int nrow, const int ncol, const int idx,
			     const int *excludes, const int m)
{
	long double rowmax = __impf_NINF__;
	int j;

	for (j = 0; j < ncol; j++) {
		if (is_in_arri(j, excludes, m))
			continue;
		if (__impf_ABS__(table[j + idx * ldtable]) > rowmax)
			rowmax = __impf_ABS__(table[j + idx * ldtable]);
	}
	return rowmax;
}

long double table_abs_colmax(const long double *table, const int ldtable,
			     const int nrow, const int ncol, const int idx)
{
	long double colmax = __impf_NINF__;
	int i;

	for (i = 0; i < nrow; i++) {
		if (__impf_ABS__(table[idx + i * ldtable]) > colmax)
			colmax = __impf_ABS__(table[idx + i * ldtable]);
	}
	return colmax;
}

/* Find the minimum non-zero absolute value of a row
 *
 * Note: a zero identifier `zeroidf` is required
 */
long double table_abs_rowmin(const long double *table, const int ldtable,
			     const int nrow, const int ncol, const int idx,
			     const int *excludes, const int m,
			     const long double zeroidf)
{
	long double e, rowmin = __impf_INF__;
	int j;

	for (j = 0; j < ncol; j++) {
		if (is_in_arri(j, excludes, m))
			continue;
		e = __impf_ABS__(table[j + idx * ldtable]);
		if (e > zeroidf && e < rowmin)
			rowmin = e;
	}
	return rowmin;
}
