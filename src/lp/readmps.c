/*
 * Copyright (C) 2022 Zhuang Linsheng <zhuanglinsheng@outlook.com>
 * License: LGPL 3.0 <https://www.gnu.org/licenses/lgpl-3.0.html>
 */
#include <impf/fmin_lp.h>
#include <stdio.h>

static struct impf_Model_LP *create_model(const int m, const int n)
{
	double *obj = NULL;
	double *coefficients = NULL;
	struct impf_LinearConstraint *constraints = NULL;
	struct impf_VariableBound *bounds = NULL;
	struct impf_Model_LP *model = NULL;
	int i;

	model = impf_malloc(sizeof(struct impf_Model_LP));
	if (model == NULL)
		return NULL;
	obj = impf_malloc(n * sizeof(double));
	if (obj == NULL) {
		impf_free(model);
		return NULL;
	}
	coefficients = impf_malloc(m * n * sizeof(double));
	if (coefficients == NULL) {
		impf_free(model);
		impf_free(obj);
		return NULL;
	}
	constraints = impf_malloc(m * sizeof(struct impf_LinearConstraint));
	if (constraints == NULL) {
		impf_free(model);
		impf_free(obj);
		impf_free(coefficients);
		return NULL;
	}
	bounds = impf_malloc(n * sizeof(struct impf_VariableBound));
	if (bounds == NULL) {
		impf_free(model);
		impf_free(obj);
		impf_free(coefficients);
		impf_free(constraints);
		return NULL;
	}
	model->m = m;
	model->n = n;
	model->objective = obj;
	model->coefficients = coefficients;
	model->constraints = constraints;
	model->bounds = bounds;
	impf_memset(coefficients, 0., m * n * sizeof(double));

	for (i = 0; i < m; i++) {
		constraints[i].coef = coefficients + i * n;
		impf_memset(constraints[i].name, '\0', 16);
	}
	for (i = 0; i < n; i++) {
		impf_memset(bounds[i].name, '\0', 16);
		bounds[i].lb = 0;
		bounds[i].ub = __impf_INF__;
		bounds[i].b_type = impf_BOUND_T_LO;
		bounds[i].v_type = impf_VAR_T_REAL;
	}
	return model;
}

void impf_lp_free(struct impf_Model_LP *model)
{
	if (model == NULL)
		return;
	if (model->coefficients)
		impf_free(model->coefficients);
	if (model->constraints)
		impf_free(model->constraints);
	if (model->bounds)
		impf_free(model->bounds);
	if (model->objective)
		impf_free(model->objective);
	impf_free(model);
}

void file_readline(FILE *f, char *line, const int n)
{
	impf_memset(line, '\0', n);
	fgets(line, n, f);
	line[impf_strcspn(line, "\r\n")] = '\0';
}

static int change_sect_code(const char *line, int *sect_code)
{
	int old_code = *sect_code;

	if (impf_memcmp(line, "ROWS", 4) == 0)
		*sect_code = 1;
	if (impf_memcmp(line, "COLUMNS", 7) == 0)
		*sect_code = 2;
	if (impf_memcmp(line, "RHS", 3) == 0)
		*sect_code = 3;
	if (impf_memcmp(line, "RANGES", 6) == 0)
		*sect_code = 4;
	if (impf_memcmp(line, "BOUNDS", 6) == 0)
		*sect_code = 5;
	if (impf_memcmp(line, "ENDATA", 6) == 0)
		*sect_code = 9;
	if (old_code == *sect_code)
		return 0;
	else
		return 1;
}

static int get_mps_info(const char *file, int *nrow, int *nvar)
{
	char line[128];
	char last_var[8];
	int sect_code = 0;

	FILE *f = fopen(file, "r");

	if (f == NULL) {
		printf("Cannot open file: \"%s\"\n", file);
		return impf_EXIT_FAILURE;
	}
	*nrow = 0;
	*nvar = 0;
	impf_memset(last_var, '\0', 8);
LOOP:
	file_readline(f, line, 128);

	if (change_sect_code(line, &sect_code))
		goto LOOP;
	switch (sect_code) {
	case 1:
		(*nrow)++;
		break;
	case 2:
		if (impf_memcmp(last_var, line + 4, 8) != 0) {
			(*nvar)++;
			impf_memcpy(last_var, line + 4, 8);
		}
		break;
	default:
		break;
	}
	if (feof(f))
		goto END;
	goto LOOP;
END:
	fclose(f);
	return impf_EXIT_SUCCESS;
}

static double get_filed_1_value(const char *line)
{
	char value_str[13];
	double value;

	impf_memset(value_str, '\0', 12);
	impf_memcpy(value_str, line + 24, 12);
	value = impf_atof(value_str);
	return value;
}

static double get_field_2_value(const char *line)
{
	char value_str[13];
	double value;

	impf_memset(value_str, '\0', 12);
	impf_memcpy(value_str, line + 49, 12);
	value = impf_atof(value_str);
	return value;
}

static void fill_model_coef(struct impf_Model_LP *model, const double value, const char *field_name, const int nvars)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (impf_memcmp(field_name, tmp, impf_strlen(tmp)) == 0) {
			model->constraints[i].coef[nvars - 1] = value;
			break;
		}
	}
}

static void fill_columns_to_model(struct impf_Model_LP *model, const char *obj_name, const char *field_name,
				  const double value, const int nvars)
{
	if (impf_memcmp(field_name, obj_name, impf_strlen(obj_name)) == 0)
		model->objective[nvars - 1] = value;
	else
		fill_model_coef(model, value, field_name, nvars);
}

static void fill_model_rhs(struct impf_Model_LP *model, const char *field_name, const double value)
{
	int i, m = model->m;

	for (i = 0; i < m; i++) {
		char *tmp = model->constraints[i].name;

		if (impf_memcmp(field_name, tmp, impf_strlen(tmp)) == 0) {
			model->constraints[i].rhs = value;
			break;
		}
	}
}

static int fill_model(const char *file, struct impf_Model_LP *model)
{
	char line[128];
	char obj_name[9];
	char *last_name = model->bounds->name;
	int sect_code = 0;
	int ncons = 0, nvars = 0;
	double value;

	FILE *f = fopen(file, "r");

	if (f == NULL) {
		printf("Cannot open file: \"%s\"\n", file);
		return impf_EXIT_FAILURE;
	}
	impf_memset(line, '\0', 128);
	impf_memset(obj_name, '\0', 9);
	impf_memset(last_name, '\0', 16);
LOOP:
	file_readline(f, line, 128);

	if (change_sect_code(line, &sect_code))
		goto LOOP;
	switch (sect_code) {
	case 1:  /* ROWS */
		switch (line[1]) {
		case 'N':
			impf_memset(obj_name, '\0', 8);
			impf_memcpy(obj_name, line + 4, 8);
			break;
		case 'L':
			impf_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = impf_CONS_T_LE;
			ncons++;
			break;
		case 'G':
			impf_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = impf_CONS_T_GE;
			ncons++;
			break;
		case 'E':
			impf_memcpy(model->constraints[ncons].name, line + 4, 8);
			model->constraints[ncons].type = impf_CONS_T_EQ;
			ncons++;
			break;
		default:
			break;
		}
		break;
	case 2:  /* COLUMNS */
		if (impf_memcmp(last_name, line + 4, 8) != 0) {
			impf_memset(model->bounds[nvars].name, '\0', 16);
			impf_memcpy(model->bounds[nvars].name, line + 4, 8);
			last_name = model->bounds[nvars].name;
			nvars++;
		}
		value = get_filed_1_value(line);
		fill_columns_to_model(model, obj_name, line + 14, value, nvars);
		if (impf_strlen(line) < 40)
			goto LOOP;
		value = get_field_2_value(line);
		fill_columns_to_model(model, obj_name, line + 39, value, nvars);
		break;
	case 3:  /* RHS */
		value = get_filed_1_value(line);
		fill_model_rhs(model, line + 14, value);
		if (impf_strlen(line) < 40)
			goto LOOP;
		value = get_field_2_value(line);
		fill_model_rhs(model, line + 39, value);
		break;
	default:
		break;
	}
	if (feof(f))
		goto END;
	goto LOOP;
END:
	fclose(f);
	return impf_EXIT_SUCCESS;
}

struct impf_Model_LP *impf_lp_readmps(const char *file)
{
	struct impf_Model_LP *model;
	int m, n;  /* number of constraints and variables */
	int n_sect_row = 0, n_sect_columns = 0;

	if (get_mps_info(file, &n_sect_row, &n_sect_columns) == impf_EXIT_FAILURE)
		return NULL;
	m = n_sect_row - 1;  /* the objective is also counted */
	n = n_sect_columns;
	model = create_model(m, n);

	if (model == NULL)
		return NULL;
	fill_model(file, model);
	return model;
}
