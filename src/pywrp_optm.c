#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#include <Python.h>
#include <impf/fmin.h>
#include <impf/fmin_lp.h>

/*
 * (x, code) = wrapper_impf_lp_simplex(m, n, maxiter, method, bounds, obj, consts...)
 */
static PyObject* wrapper_impf_lp_simplex(PyObject* self, PyObject* args) {
	int m, n, max_iter, i, j;
	char *method;
	PyObject *bounds = NULL; /* python list of tuple */
	PyObject *obj; /* python list of float */
	PyObject *consts_coef; /* python list of list of float */
	PyObject *consts_rhs; /* python list of float */
	PyObject *consts_type; /* python list of int */

	if(!PyArg_ParseTuple(args, "iiis|OOOOO",&m, &n, &max_iter, &method,
		&bounds, &obj, &consts_coef, &consts_rhs, &consts_type)) {
		return NULL;
	}

	/* After assign m and n */
	double *obj_val = impf_malloc(sizeof(double) * n);
	struct impf_VariableBound *bounds_val = NULL;
	struct impf_LinearConstraint *constraints = impf_malloc(sizeof(struct impf_LinearConstraint) * m);
	double *x = impf_malloc(sizeof(double) * n);
	double value;
	int error_code;

	/* Parse bounds */
	if (bounds != Py_None) {
		bounds_val = impf_malloc(sizeof(struct impf_VariableBound) * n);

		for (j = 0; j<n; j++) {
			PyObject *item = PyList_GetItem(bounds, j);
			PyObject *lb_obj = PyTuple_GetItem(item, 0);
			PyObject *ub_obj = PyTuple_GetItem(item, 1);
			double lb = PyFloat_AsDouble(lb_obj);
			double ub = PyFloat_AsDouble(ub_obj);

			if (1.0 / lb == 0.0 && lb < 0)
				bounds_val[j].lb = __impf_NINF__;
			if (1.0 / ub == 0.0 && ub > 0)
				bounds_val[j].ub = __impf_INF__;
			if (bounds_val[j].ub == __impf_INF__ && bounds_val[j].lb == __impf_NINF__)
				bounds_val[j].b_type = impf_BOUND_T_FR;
			else if (bounds_val[j].ub == __impf_INF__ && bounds_val[j].lb != __impf_NINF__)
				bounds_val[j].b_type = impf_BOUND_T_UP;
			else if (bounds_val[j].ub != __impf_INF__ && bounds_val[j].lb == __impf_NINF__)
				bounds_val[j].b_type = impf_BOUND_T_LO;
			else if (bounds_val[j].ub != __impf_INF__ && bounds_val[j].lb != __impf_NINF__)
				bounds_val[j].b_type = impf_BOUND_T_BS;
		}
	}

	/* Parse obj */
	for(j = 0; j < n; j++) {
		PyObject *item = PyList_GetItem(obj, j);
		obj_val[j] = PyFloat_AsDouble(item);
	}

	/* Parse consts */
	for(i = 0; i < m; i++) {
		PyObject *const_i_coef = PyList_GetItem(consts_coef, i);
		PyObject *const_i_rhs = PyList_GetItem(consts_rhs, i);
		PyObject *const_i_type = PyList_GetItem(consts_type, i);
		double *const_i_coef_value = impf_malloc(sizeof(double) * n);

		for(j = 0; j < n; j++) {
			PyObject *item = PyList_GetItem(const_i_coef, j);
			const_i_coef_value[j] = PyFloat_AsDouble(item);
		}
		constraints[i].coef = const_i_coef_value;
		constraints[i].rhs = PyFloat_AsDouble(const_i_rhs);
		constraints[i].type = (int) PyLong_AsLong(const_i_type);
	}

	/* Call solver */
	impf_lp_simplex(obj_val, constraints, bounds_val, m, n, method, max_iter, x, &value, &error_code);

	/* Return object */
	PyObject *py_x_list = PyList_New(n);
	PyObject* result = PyTuple_New(2);

	for (j = 0; j < n; j++) {
		PyObject *item = PyFloat_FromDouble(x[j]);
		PyList_SetItem(py_x_list, j, item);
	}
	if (x)
		impf_free(x);
	if (constraints) {
		for (i = 0; i<m; i++) {
			if (constraints[i].coef)
				impf_free(constraints[i].coef);
		}
		impf_free(constraints);
	}
	if (obj_val)
		impf_free(obj_val);
	if (bounds_val)
		impf_free(bounds_val);
	PyTuple_SetItem(result, 0, py_x_list);
	PyTuple_SetItem(result, 1, PyLong_FromLong(error_code));
	return result;
}

/* Module method table */
static PyMethodDef optm_Methods[] = {
	{"wrapper_impf_lp_simplex", wrapper_impf_lp_simplex, METH_VARARGS, "Linear Programming"},
	{ NULL, NULL, 0, NULL}
};

/* Module structure */
static struct PyModuleDef optm_module = {
	PyModuleDef_HEAD_INIT,
	"_clib_optm",
	"Optimization Module",
	-1,
	optm_Methods
};

/* Module initialization function */
PyMODINIT_FUNC PyInit__clib_optm(void) {
	return PyModule_Create(&optm_module);
}
