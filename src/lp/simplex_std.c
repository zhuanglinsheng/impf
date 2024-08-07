/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/_magics.h>
#include <impf/linalg.h>
#include <impf/fmin_lp.h>
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
 *     0    not circled
 *     1    circled but optimal already
 *     2    circled
 */
static int check_simplex_degenerated(const double *table, const int n, const double old_value)
{
	if (old_value <= table[n] + __impf_CHC_SPLX_DEGENERATED__) {
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

		if (y_i_q <= __impf_CTR_SPLX_PIV_LEV__)
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
static int simplex_pivot_enter_rule_datzig(const double *table, const int *basis, const int m, const int n)
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
static int simplex_pivot_enter_rule_bland(const double *table, const int *basis, const int m, const int n)
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
			       const int rule1, const int rule2, const int rule3)
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
 *	0: current BFS is NOT optimal (stop before converged)
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
	if (7 == impf_strlen(criteria) && 0 == impf_memcmp("dantzig", criteria, 7)) {
		q = simplex_pivot_enter_rule_datzig(table, basis, m, n);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);
	}
	else if (5 == impf_strlen(criteria) && 0 == impf_memcmp("bland", criteria, 5)) {
		q = simplex_pivot_enter_rule_bland(table, basis, m, n);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);
	}
	else {  /* default method: "pan97" */
		q = simplex_pivot_enter_rule_datzig(table, basis, m, n);
		p = simplex_pivot_leave_rule(table, ldtable, m, n, q, &bounded);

		/* int idx_degen = find_bv_degenerated(table, ldtable, m, n); */
#ifdef IMPF_MODE_DEBUG
		impf_prt_matd(table, n + 1, m + 1, n + 1);
		printf("\n");
		printf("p = %i, q = %i\nbasis = ", p, q);
		impf_prt_arri(basis, m, 1);
		printf("\n");
#endif
	}
	if (n <= q) {
#ifdef IMPF_MODE_DEBUG
		printf("Pivot failure due to '9: numerical precision error'\n");
#endif
		return 9;
	}

	if (bounded == 0)
		return 2;
	basis[p] = q;
	simplex_pivot_core(table, ldtable, m, n, p, q, 1, 1, 1);
	return 0;
}

/* Linear Programming: simplex algorithm for solving LP of basic representation
 *
 * Return
 *	0: current BFS is NOT optimal (stop before converged)
 *	1: current BSF is optimal
 *	2: LP is unbounded
 *	3: LP is circled more than accepted times (indicating for degeneracy)
 *	9: numerical precision error
 */
static int simplex_pivot_bsc(int *epoch, double *table, const int ldtable, int *basis,
			     const int m, const int n, const int nreal,
			     const char *criteria, const int niter)
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
#ifdef IMPF_MODE_DEBUG
		printf(">>> Algorithm stop due to '1: current BSF is optimal'.\n");
#endif
			return 1;
		case 2:
#ifdef IMPF_MODE_DEBUG
		printf(">>> Algorithm stop due to '2: LP is unbounded'.\n");
#endif
			return 2;
		case 9:
#ifdef IMPF_MODE_DEBUG
		printf(">>> Algorithm stop due to '9: numerical precision error'.\n");
#endif
			return 9;
		}
		if (check_simplex_degenerated(table, n, old_value) == 2) {
#ifdef IMPF_MODE_DEBUG
			printf(">>> Degenerated, value = %f [%i]\n", table[n], degen_iter);
#endif
			degen_iter++;
			if (degen_iter > 5)
				return 3;
		} else
			degen_iter = 0;
		old_value = table[n];
	}
	return 0;
}

/* To create in heap (need to be released) simplex table, index set of basis
 * and constraint type recorder
 */
static int simplex_create_buffer(double **table, int **basis, int **constypes,
				 const int m, const int nrow, const int ncol)
{
	*table = NULL;
	*basis = NULL;
	*constypes = NULL;

	*table = impf_malloc(nrow * ncol * sizeof(double));
	if (*table == NULL)
		return impf_EXIT_FAILURE;
	*basis = impf_malloc(m * sizeof(int));
	if (*basis == NULL) {
		impf_free(*table);
		return impf_EXIT_FAILURE;
	}
	*constypes = impf_malloc(m * sizeof(int));
	if (*constypes == NULL) {
		impf_free(*table);
		impf_free(*basis);
		return impf_EXIT_FAILURE;
	}
	return impf_EXIT_SUCCESS;
}

static void simplex_free_buffer(double *table, int *basis, int *constypes)
{
	if (table)
		impf_free(table);
	if (basis)
		impf_free(basis);
	if (constypes)
		impf_free(constypes);
}

/* Fill in constraint type array from "constraints"
 *
 * Constraints rhs are transformed to be nonnegative,
 * "LE" and "GE" types are transformed respectively
 */
static void simplex_fill_constypes(const struct impf_LinearConstraint *constraints, int *constypes, const int m)
{
	int i;

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;

		if (cons->rhs >= 0)
			constypes[i] = cons->type;
		else {
			switch (cons->type) {
			case impf_CONS_T_EQ:
				constypes[i] = impf_CONS_T_EQ;
				break;
			case impf_CONS_T_GE:
				constypes[i] = impf_CONS_T_LE;
				break;
			case impf_CONS_T_LE:
				constypes[i] = impf_CONS_T_GE;
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

	impf_memset(table, 0., nrow * ldtable * sizeof(double));

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;
		int row = (i + 1) * ldtable;

		if (cons->rhs >= 0) {
			table[ncol - 1 + row] = cons->rhs;
			impf_memcpy(table + row, cons->coef, n * sizeof(double));
		} else {
			table[ncol - 1 + row] = -cons->rhs;
			for (j = 0; j < n; j++)
				table[j + row] = -cons->coef[j];
		}
	}
}

/* Determine the size of (basic) simplex table
 *
 * "GE" constraint has a slack var and an artificial var, hence will generate
 * an additional variable than usual
 */
static void simplex_table_size_usul(const struct impf_LinearConstraint *constraints,
				    const int m, const int n, int *nrow, int *ncol)
{
	int i;
	*nrow = m + 1;
	*ncol = m + n + 1;

	for (i = 0; i < m; i++) {
		const struct impf_LinearConstraint *cons = constraints + i;

		if (impf_CONS_T_GE == cons->type && cons->rhs >= 0)
			(*ncol)++;
		if (impf_CONS_T_LE == cons->type && cons->rhs < 0)
			(*ncol)++;
	}
}

/* Add slack variables (GE, LE) to simplex table
 * Return the number of slack variables
 */
static int simplex_add_slack(double *table, const int ldtable, const int *constypes, const int m, const int n)
{
	int i, nslack = 0;

	for (i = 0; i < m; i++) {
		if (impf_CONS_T_GE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = -1.;
			nslack++;
		}
		if (impf_CONS_T_LE == constypes[i]) {
			table[n + nslack + (i + 1) * ldtable] = 1.;
			nslack++;
		}
	}
	return nslack;
}

/* Add artificial variables (GE, EQ) to simplex table
 * Return the number of artificial variables
 */
static int simplex_add_artif(double *table, const int ldtable, const int *constypes,
			     const int m, const int n, const int nslack)
{
	int i, nartif = 0;

	for (i = 0; i < m; i++) {
		if (impf_CONS_T_LE != constypes[i]) {
			table[n + nslack + nartif + (i + 1) * ldtable] =  1.;
			table[n + nslack + nartif] = -1.;
			nartif++;
		}
	}
	return nartif;
}

/* Fill in the basis index set of artificial LP
 */
static void simplex_fill_artiflp_basis(int *basis, const int *constypes,
				       const int m, const int n, const int nslack)
{
	int i, tmp_nbasis = 0, tmp_nslack = 0, tmp_nartif = 0;

	for (i = 0; i < m; i++) {
		switch (constypes[i]) {
		case impf_CONS_T_EQ:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			break;
		case impf_CONS_T_GE:
			*(basis + tmp_nbasis) = n + nslack + tmp_nartif;
			tmp_nartif++;
			tmp_nslack++;
			break;
		case impf_CONS_T_LE:
			*(basis + tmp_nbasis) = n + tmp_nslack;
			tmp_nslack++;
			break;
		}
		tmp_nbasis++;
	}
}

static void simplex_fill_artiflp_nrcost(double *table, const int ldtable, const int *constypes,
					const int m, const int ncol)
{
	int i, rowi;

	for (i = 0; i < m; i++) {
		if (impf_CONS_T_LE == constypes[i])
			continue;
		rowi = (i + 1) * ldtable;
		impf_linalg_daxpy(ncol, 1, table + rowi, 1, table, 1);
	}
}

static void transf_artif_basis(double *table, int ldtable, int *basis, const int m, const int nreal, int nvar)
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
			impf_memset(table + (i + 1) * ldtable, 0, nreal * sizeof(double));
			continue;
		}
		simplex_pivot_core(table, ldtable, m, nvar, i, q, 1, 1, 0);
		basis[i] = q;
#ifdef IMPF_MODE_DEV
		printf("in 'transf_artif_basis'\n");
		printf("switch %i <----> %i\n", i, q);
#endif
		q = nvar;
	}
}

static void delete_artif_cols(double *table, const int ldtable, const int m, const int nreal, const int nartif)
{
	int i, rowi;

	if (nartif <= 0)  /* Delete artificial columns */
		return;
	for (i = 0; i < m + 1; i++) {
		rowi = i * ldtable;
		table[nreal + rowi] = table[nreal + nartif + rowi];
	}
}

/* Phase 1: get a BFS for the original problem using the usual way - artificial LP
 *
 * Work:
 * 	1. allocate memory for table, basis, constypes
 * 	2. form a basic feasible solution (BSF)
 * 	3. assign ldtable and nvar, the number of vars in BSF
 */
static int simplex_phase_1_usul(double **table, int *ldtable, int **basis, int **constypes,
				int *nvar, int *epoch, int *code,
				const struct impf_LinearConstraint *constraints,
				const int m, const int n, const char *criteria, const int niter)
{
	int nrow, ncol;
	int nslack, nartif;

#ifdef IMPF_MODE_DEV
	printf(">>> Reformulate, code = %i\n", *code);
#endif
	simplex_table_size_usul(constraints, m, n, &nrow, &ncol);
	*ldtable = ncol;  /* leading dimension of table in memory */
	if (simplex_create_buffer(table, basis, constypes, m, nrow, *ldtable) == impf_EXIT_FAILURE) {
		*code = impf_MemoryAllocError;
		return impf_EXIT_FAILURE;
	}

	simplex_fill_constypes(constraints, *constypes, m);
	simplex_fill_conscoefs(*table, *ldtable, constraints, nrow, ncol, m, n);
	nslack = simplex_add_slack(*table, *ldtable, *constypes, m, n);
	nartif = simplex_add_artif(*table, *ldtable, *constypes, m, n, nslack);
	*nvar = n + nslack + nartif;  /* will be recovered to `n + nslack` upon success */
	if (m > (*nvar)) {
		*code = impf_OverDetermination;
		goto END;
	}
	simplex_fill_artiflp_basis(*basis, *constypes, m, n, nslack);
	simplex_fill_artiflp_nrcost(*table, *ldtable, *constypes, m, ncol);

#ifdef IMPF_MODE_DEV
	printf(">>> n = %i, nslack = %i, nartif = %i\n", n, nslack, nartif);
	printf(">>> table size = (%i, %i)\n", nrow, ncol);
	printf(">>> Pivoting, code = %i\n", *code);
#endif
	switch (simplex_pivot_bsc(epoch, *table, *ldtable, *basis, m, *nvar, n + nslack, criteria, niter)) {
	case 0:
		*code = impf_ExceedIterLimit;
		goto END;
	case 1:
		if ((*table)[ncol - 1] > __impf_CHC_SPLX_FEASIBLE__) {
			*code = impf_Infeasibility;
			goto END;
		}
		transf_artif_basis(*table, *ldtable, *basis, m, n + nslack, *nvar);
		delete_artif_cols(*table, *ldtable, m, n + nslack, nartif);
		*nvar = n + nslack;
		return impf_EXIT_SUCCESS;
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
	simplex_free_buffer(*table, *basis, *constypes);
	return impf_EXIT_FAILURE;
}

/* Phase 2: Solve the original problem
 */
static int simplex_phase_2_usul(double *table, int ldtable, int *basis, int *constypes,
				int *epoch, int *code, const int m, const int n,
				const int nvar, const char *criteria, const int niter)
{
#ifdef IMPF_MODE_DEV
	printf(">>> m = %i, n = %i, nvar = %i\n", m, n, nvar);
	printf(">>> Pivoting, code = %i\n", *code);
#endif
	switch (simplex_pivot_bsc(epoch, table, ldtable, basis, m, nvar, nvar, criteria, niter)) {
	case 0:
		*code = impf_ExceedIterLimit;
		goto END;
	case 1:
		*code = impf_Success;
		return impf_EXIT_SUCCESS;
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
	return impf_EXIT_FAILURE;  /* error code already updated */
}

int impf_lp_simplex_std(const double *objective, const struct impf_LinearConstraint *constraints,
			const int m, const int n, const char *criteria, const int niter,
			double *x, double *value, int *code)
{
	int i, j;
	int ldtable;
	int nvar;
	int epoch = 0;
	int *basis = NULL;
	double *table = NULL;
	int *constypes = NULL;

	assert(objective != NULL);
	assert(constraints != NULL);
	assert(x != NULL);
	assert(value != NULL);
	assert(code != NULL);

#ifdef IMPF_MODE_DEV
	printf("Phase 1 Begin: code = %i\n", *code);
#endif
	if (simplex_phase_1_usul(&table, &ldtable, &basis, &constypes, &nvar, &epoch, code,
				 constraints, m, n, criteria, niter) == impf_EXIT_FAILURE)
		return impf_EXIT_FAILURE;
#ifdef IMPF_MODE_DEV
	printf("Phase 1 Done.\n");
#endif

#ifdef IMPF_MODE_DEV
	printf("Phase 2 Begin:\n");
#endif
	for (j = 0; j < n; j++)  /* Fill in original objective coefficients */
		table[j] = -objective[j];
	for (i = 0; i < m; i++) {  /* row_0 = row_0 - ratio * row_{i+1} */
		int rowi = (i + 1) * ldtable;
		double ratio = -table[basis[i]];

		impf_linalg_daxpy(nvar + 1, ratio, table + rowi, 1, table, 1);
	}
	if (simplex_phase_2_usul(table, ldtable, basis, constypes, &epoch, code,
				 m, n, nvar, criteria, niter) == impf_EXIT_FAILURE)
		return impf_EXIT_FAILURE;
#ifdef IMPF_MODE_DEV
	printf("Phase 2 Done.\n");
#endif

	*value = table[nvar];
	impf_memset(x, 0., n * sizeof(double));
	for (i = 0; i < m; i++) {
		if (basis[i] < n)
			x[basis[i]] = table[nvar + (i + 1) * ldtable];
	}
	simplex_free_buffer(table, basis, constypes);
	return impf_EXIT_SUCCESS;
}
