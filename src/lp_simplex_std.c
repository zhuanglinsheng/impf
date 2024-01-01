/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>
#include <impf/lp.h>
#include <string.h>
#ifdef IMPF_MODE_DEV
#include <stdio.h>
#endif

static int is_simplex_optimal(const double *table, const int n)
{
	int j;

	for (j = 0; j < n; j++) {
		if (table[j] > __impf_CTR_SPLX_OPTIMAL__)
			return 0;
	}
	return 1;
}

/* Basis index should not exceeds the number of variables
 */
static int is_basis_valid(const int *basis, const int m, const int n)
{
	return maxabs_arri(basis, m, 1) <= n;
}

/* In simplex iteration, check whether value is improved.
 *
 * Return:
 *     0    not degenerated
 *     1    degenerated but optimal already
 *     2    degenerated
 */
static int check_simplex_degeneracy(const double *table, const int n, const double old_value)
{
	if (old_value <= table[n] + __impf_CHC_SPLX_DEGEN__) {
		if (is_simplex_optimal(table, n))
			return 1;
		return 2;
	} else
		return 0;
}

/* Choose the variable to leave basis
 * Return the index of the variable and check weather LP is "bounded"
 */
static int simplex_pivot_leave_rule(const double *table, const int ldtable,
				    const int m, const int n, const int q, int *bounded)
{
	int i, p = n;
	double y_i_0, y_i_q, x_iq, min_x_iq = __impf_INF__;
	*bounded = 0;

	for (i = 0; i < m; i++) {
		y_i_0 = table[n + (i + 1) * ldtable];
		y_i_q = table[q + (i + 1) * ldtable];

		if (y_i_q <= __impf_CTR_SPLX_PIVLEV_ZERO__)
			continue;
		else {
			x_iq = y_i_0 / y_i_q;

			if (x_iq < min_x_iq) {
				min_x_iq = x_iq;
				p = i;
			}
			*bounded = 1;
		}
	}
	return p;
}

/* Fast pivot rule: choosing the variable to enter basis
 * Return the index of the variable (< n)
 *
 * Note:
 * 	On failure, the algorithm returns `n`. Logically, this NEVER happens,
 *	but numerically, there are many criteriors reporting optimality,
 *	leading to unpredicted results
 */
static int simplex_dantzig_enter_rule(const double *table, const int *basis, const int m, const int n)
{
	int j, q = n;
	double beta_j = 0.;
	double beta_q = 0.;

	for (j = 0; j < n; j++) {
		if (is_in_arri(j, basis, m))
			continue;
		beta_j = table[j];

		if (beta_j > beta_q) {
			q = j;
			beta_q = beta_j;
		}
	}
	return q; /* return n if no p is found (optimal already) */
}

/* Bland's rule: choosing the variable to enter basis
 * Return the index of the variable (< n)
 *
 * Note: on failure, the algorithm returns n
 */
static int simplex_bland_enter_rule(const double *table, const int *basis, const int m, const int n)
{
	int j;
	double epsilon = __impf_CTR_SPLX_BLAND_EPS__;
BLAND_BEGIN:
	for (j = 0; j < n; j++) {
		if (!is_in_arri(j, basis, m)) {
			if (table[j] > epsilon)
				return j;
		}
	}
	if (epsilon >= __impf_CTR_SPLX_BLAND_EPS_MIN__) {
		epsilon /= 10.;
		goto BLAND_BEGIN;
	} else
		return n;
}

/* Key subroutine of pivoting
 *
 * Parameter:
 *	p	idx of variable to leave basis
 *	q	idx of variable to enter basis
 *
 * Work:
 *	rule 1. row_p normalized by dividing y_p_q
 *	rule 2. row_i -= row_p * y_i_q
 *	rule 3. row_0 -= row_p * beta_q
 */
static void simplex_pivot_core(double *table, const int ldtable,
			       const int m, const int n, const int p, const int q,
			       int rule1, int rule2, int rule3)
{
	int i, ncol = n + 1, rowp = (p + 1) * ldtable;
	double y_p_q = table[q + rowp];

	if (rule1)
		impf_linalg_dscal(ncol, 1 / y_p_q, table + rowp, 1);
	if (rule2) {
		for (i = 0; i < m; i++) {
			int rowi = (i + 1) * ldtable;
			double rto =  -table[q + rowi];

			if (i == p)
				continue;
			impf_linalg_daxpy(ncol, rto, table + rowp, 1, table + rowi, 1);
		}
	}
	if (rule3)
		impf_linalg_daxpy(ncol, -table[q], table + rowp, 1, table, 1);
}

/* Pivot starting from a basic representation for one round
 *
 * Return:
 *	0: current BFS is NOT optimal
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	9: numerical precision error
 */
static int simplex_pivot_on(double *table, const int ldtable, int *basis,
			    const int m, const int n, const char *criteria)
{
	int bounded;
	int q, p;

	if (is_simplex_optimal(table, n))
		return 1;
	if (7 == strlen(criteria) && 0 == memcmp("dantzig", criteria, 7))
		q = simplex_dantzig_enter_rule(table, basis, m, n);
	else if (5 == strlen(criteria) && 0 == memcmp("bland", criteria, 5))
		q = simplex_bland_enter_rule(table, basis, m, n);
	else {  /* default method: "pan97" */
		return 9;
	}
	if (n <= q)
		return 9;
	p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);

	if (bounded == 0)
		return 2;
	basis[p] = q;
	simplex_pivot_core(table, ldtable, m, n, p, q, 1, 1, 1);
	return 0;
}

/* Linear Programming: simplex algorithm for solving LP of basic representation
 *
 * Return
 *	0: current BFS is NOT optimal
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	3: LP is degenerated
 *	9: numerical precision error
 */
static int simplex_pivot_bsc(int *epoch, double *table, const int ldtable, int *basis,
			     const int m, const int n, const int nreal,
			     const char *criteria, const int niter, const int phase1)
{
	double old_value = __impf_INF__;
	int degen_iter = 0;

	assert(table != NULL);
	assert(basis != NULL);
	assert(epoch != NULL);

	while (*epoch < niter) {
		(*epoch)++;
		switch (simplex_pivot_on(table, ldtable, basis, m, n, criteria)) {
		case 0:
			break;
		case 1:
			return 1;
		case 2:
			return 2;
		case 9:
			return 9;
		}
#ifdef IMPF_MODE_DEBUG
		printf("[%u] value = %f - ", *epoch, table[n]);
		printf("degen-iter = %u\n\n", degen_iter);
#endif
		if (check_simplex_degeneracy(table, n, old_value) == 2) {
			degen_iter++;
			if (degen_iter > 50)
				return 3;
		} else
			degen_iter = 0;
		old_value = table[n];

		if (phase1 && table[n] < -__impf_CHC_SPLX_PHASE_1_NNVAL__) {
#ifdef IMPF_MODE_DEV
			printf("NEGATIVE VALUE IN PHASE 1\n");
#ifdef IMPF_MODE_DEBUG
			printf("[%u] value = %f\n", *epoch, table[n]);
			printf("optimal = %u\n", is_simplex_optimal(table, n));
#endif
#endif
			return 9;
		}
	}
	return 0;
}

/* Determine the size of (basic) simplex table
 *
 * "GE" constraint has a slack var and an artificial var, hence will generate
 * an additional variable than usual
 */
static void simplex_table_size(const struct impf_LinearConstraint *constraints,
			       const int m, const int n, int *nrow, int *ncol)
{
	int i;
	*nrow = m + 1;
	*ncol = m + n + 1;

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;

		if (impf_GE == cons->type && cons->rhs >= 0)
			(*ncol)++;
		if (impf_LE == cons->type && cons->rhs < 0)
			(*ncol)++;
	}
}

/* To create in heap (need to be released) simplex table, index set of basis
 * and constraint type recorder
 */
static int simplex_create_buffer(double **table, int **basis, enum impf_ConstraintType **constypes,
				 const int m, const int nrow, const int ncol)
{
	*table = NULL;
	*basis = NULL;
	*constypes = NULL;

	*table = malloc(nrow * ncol * sizeof(double));
	if (*table == NULL)
		return EXIT_FAILURE;
	*basis = malloc(m * sizeof(int));
	if (*basis == NULL) {
		free(*table);
		return EXIT_FAILURE;
	}
	*constypes = malloc(m * sizeof(enum impf_ConstraintType));
	if (*constypes == NULL) {
		free(*table);
		free(*basis);
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}

static void simplex_free_buffer(double *table, int *basis, enum impf_ConstraintType *constypes)
{
	if (table)
		free(table);
	if (basis)
		free(basis);
	if (constypes)
		free(constypes);
}

/* Fill in constraint type array from "constraints"
 *
 * Constraints rhs are transformed to be nonnegative,
 * "LE" and "GE" types are transformed respectively
 */
static void simplex_fill_constypes(const struct impf_LinearConstraint *constraints,
				   enum impf_ConstraintType *constypes, const int m)
{
	int i;

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;

		if (cons->rhs >= 0)
			constypes[i] = cons->type;
		else {
			switch (cons->type) {
			case impf_EQ:
				constypes[i] = impf_EQ;
				break;
			case impf_GE:
				constypes[i] = impf_LE;
				break;
			case impf_LE:
				constypes[i] = impf_GE;
				break;
			}
		}
	}
}

/* Fill in coef and rhs of constraints
 *
 * Constraints rhs are transformed to be nonnegative
 */
static void simplex_fill_conscoefs(double *table, const int ldtable, const struct impf_LinearConstraint *constraints,
				   const int nrow, const int ncol, const int m, const int n)
{
	int i, j;

	memset(table, 0., nrow * ldtable * sizeof(double));

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;
		int row = (i + 1) * ldtable;

		if (cons->rhs >= 0) {
			table[ncol - 1 + row] = cons->rhs;
			memcpy(table + row, cons->coef, n * sizeof(double));
		} else {
			table[ncol - 1 + row] = -cons->rhs;
			for (j = 0; j < n; j++)
				table[j + row] = -cons->coef[j];
		}
	}
}

/* Add slack variables (GE, LE) to simplex table
 * Return the number of slack variables
 */
static int simplex_add_slack(double *table, const int ldtable, const enum impf_ConstraintType *constypes,
			     const int m, const int n)
{
	int i, nslack = 0;

	for (i = 0; i < m; i++) {
		if (impf_GE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = -1.;
			nslack++;
		}
		if (impf_LE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = 1.;
			nslack++;
		}
	}
	return nslack;
}

/* Add artificial variables (GE, EQ) to simplex table
 * Return the number of artificial variables
 */
static int simplex_add_artif(double *table, const int ldtable, const enum impf_ConstraintType *constypes,
			     const int m, const int n, const int nslack)
{
	int i, nartif = 0;

	for (i = 0; i < m; i++) {
		if (impf_LE != constypes[i]) {
			table[n + nslack + nartif + (i + 1) * ldtable] =  1.;
			table[n + nslack + nartif] = -1.;
			nartif++;
		}
	}
	return nartif;
}

/* Fill in the basis index set of artificial LP
 */
static void simplex_fill_artiflp_basis(int *basis, const enum impf_ConstraintType *constypes,
				       const int m, const int n, const int nslack)
{
	int i, tmp_nbasis = 0, tmp_nslack = 0, tmp_nartif = 0;

	for (i = 0; i < m; i++) {
		switch (constypes[i]) {
		case impf_EQ:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			break;
		case impf_GE:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			tmp_nslack++;
			break;
		case impf_LE:
			*(basis + tmp_nbasis) = n + tmp_nslack;
			tmp_nslack++;
			break;
		}
		tmp_nbasis++;
	}
}

static void doubleransf_artificial_basis(double *table, int ldtable, int *basis,
				const int m, const int nreal, int nvar)
{
	int i, j, q = nvar;
	double ele, maxv = __impf_NINF__;

	if (is_basis_valid(basis, m, nreal))
		return;
	for (i = 0; i < m; i++) {
		if (basis[i] <nreal)
			continue;
		for (j = 0; j < nreal; j++) {
			ele = __impf_ABS__(table[j + (i + 1) * ldtable]);
			if (ele > maxv) {
				maxv = ele;
				q = j;
			}
		}
		if (maxv < 1e-9) {
			table[nvar + (i + 1) * ldtable] = 0.;
			memset(table + (i + 1) * ldtable, 0, nreal * sizeof(double));
			continue;
		}
		simplex_pivot_core(table, ldtable, m, nvar, i, q, 1, 1, 0);
		basis[i] = q;
#ifdef IMPF_MODE_DEV
		printf("in 'doubleransf_artificial_basis'\n");
		printf("switch %i <----> %i\n", i, q);
#endif
		q = nvar;
	}
}

static void simplex_fill_artiflp_nrcost(double *table, const int ldtable, const enum impf_ConstraintType *constypes,
					const int m, const int ncol)
{
	int i, rowi;

	for (i = 0; i < m; i++) {
		if (impf_LE == constypes[i])
			continue;
		rowi = (i + 1) * ldtable;
		impf_linalg_daxpy(ncol, 1, table + rowi, 1, table, 1);
	}
}

int impf_lp_simplex_std(const double *objective, const struct impf_LinearConstraint *constraints,
			const int m, const int n, const char *criteria, const int niter,
			double *x, double *value, enum impf_ErrorCode *code)
{
	int i, j;
	int nrow, ncol, ldtable;  /* simplex table size */
	int nslack, nartif, nvar;
	double *table;
	int *basis;
	enum impf_ConstraintType *constypes;
	int epoch = 0;

	assert(objective != NULL);
	assert(constraints != NULL);
	assert(x != NULL);
	assert(value != NULL);
	assert(code != NULL);

	simplex_table_size(constraints, m, n, &nrow, &ncol);
	ldtable = ncol;  /* leading dimension of table in memory */
	if (simplex_create_buffer(&table, &basis, &constypes, m, nrow, ldtable) == EXIT_FAILURE) {
		*code = impf_MemoryAllocError;
		return EXIT_FAILURE;
	}
	simplex_fill_constypes(constraints, constypes, m);
	simplex_fill_conscoefs(table, ldtable, constraints, nrow, ncol, m, n);

	/*
	 * Phase 1: Solve the artificial problem
	 */
	nslack = simplex_add_slack(table, ldtable, constypes, m, n);
	nartif = simplex_add_artif(table, ldtable, constypes, m, n, nslack);
	nvar = n + nslack + nartif;
	if (m > nvar) {
		*code = impf_OverDetermination;
		goto END;
	}
	simplex_fill_artiflp_basis(basis, constypes, m, n, nslack);
	simplex_fill_artiflp_nrcost(table, ldtable, constypes, m, ncol);
#ifdef IMPF_MODE_DEV
	printf("Phase 1:\n");
#endif
	switch (simplex_pivot_bsc(&epoch, table, ldtable, basis, m, nvar, n + nslack, criteria, niter, 1)) {
	case 0:
		*code = impf_ExceedIterLimit;
		goto END;
	case 1:
		if (table[ncol - 1] > __impf_CHC_SPLX_FEASIBLE__) {
#ifdef IMPF_MODE_DEV
			printf("value = %e\n", table[ncol - 1]);
#endif
			*code = impf_Infeasibility;
			goto END;
		}
#ifdef IMPF_MODE_DEV
		printf("Phase 1 Done.\n");
#endif
		break;
	case 2:
		*code = impf_PrecisionError;
		goto END;
	case 3:
		*code = impf_Degeneracy;
		goto END;
	case 9:
		*code = impf_PrecisionError;
		goto END;
	}

	/*
	 * Phase 2: Solve the original problem
	 */
	doubleransf_artificial_basis(table, ldtable, basis, m, n + nslack, nvar);
	for (j = 0; j < n; j++)  /* Fill in original objective coefficients */
		table[j] = -objective[j];
	if (nartif > 0) {  /* Delete artificial columns */
		for (i = 0; i < nrow; i++) {
			int rowi = i * ldtable;

			table[n + nslack + rowi] = table[ncol - 1 + rowi];
		}
	}
	nvar = n + nslack;
	ncol = nvar + 1;
	for (i = 0; i < m; i++) {  /* row_0 = row_0 - ratio * row_{i+1} */
		int rowi = (i + 1) * ldtable;
		double ratio = -table[basis[i]];

		impf_linalg_daxpy(ncol, ratio, table + rowi, 1, table, 1);
	}
#ifdef IMPF_MODE_DEV
	printf("Phase 2:\n");
#endif
	switch (simplex_pivot_bsc(&epoch, table, ldtable, basis, m, nvar, nvar, criteria, niter, 0)) {
	case 0:
		*code = impf_ExceedIterLimit;
		goto END;
	case 1:
#ifdef IMPF_MODE_DEV
		printf("Phase 2 Done.\n");
#endif
		*value = table[nvar];
		memset(x, 0., n * sizeof(double));
		for (i = 0; i < m; i++) {
			if (basis[i] < n)
				x[basis[i]] = table[nvar + (i + 1) * ldtable];
		}
		simplex_free_buffer(table, basis, constypes);
		*code = impf_Success;
		return EXIT_SUCCESS;
	case 2:
		*code = impf_Unboundedness;
		goto END;
	case 3:
		*code = impf_Degeneracy;
		goto END;
	case 9:
		*code = impf_PrecisionError;
		goto END;
	}
END:
	simplex_free_buffer(table, basis, constypes);
	return EXIT_FAILURE;  /* error code already updated */
}
