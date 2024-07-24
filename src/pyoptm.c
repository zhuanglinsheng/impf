#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#include <Python.h>
#include <impf/fmin_lp.h>

static PyObject* fmin_lp(PyObject* self, PyObject* args) {
	int m, n, max_iter, i, j;
	char* method;
	PyObject *obj_coef; /* python list of float */
	PyObject *consts_coef; /* python list of list of float */
	PyObject *consts_rhs; /* python list of float */
	PyObject *consts_type; /* python list of int */

	if(!PyArg_ParseTuple(args, "iiisOOOO",&m, &n, &max_iter, &method,
		&obj_coef, &consts_coef, &consts_rhs, &consts_type)) {
		return NULL;
	}

	/* After assign m and n */
	double *obj_coef_value = impf_malloc(sizeof(double) * n);
	struct impf_LinearConstraint *constraints = impf_malloc(sizeof(struct impf_LinearConstraint) * m);
	double *x = impf_malloc(sizeof(double) * n);
	double value;
	int error_code;

	/* Parse obj */
	for(j = 0; j < n; j++) {
		PyObject *item = PyList_GetItem(obj_coef, j);
		obj_coef_value[j] = PyFloat_AsDouble(item);
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
	impf_lp_simplex(obj_coef_value, constraints, NULL, m, n, method, max_iter, x, &value, &error_code);

	/* Return object */
	PyObject *py_x_list = PyList_New(n);

	for (j = 0; j < n; j++) {
		PyObject *item = PyFloat_FromDouble(x[j]);
		PyList_SetItem(py_x_list, j, item);
	}
	impf_free(x);
	impf_free(constraints);
	impf_free(obj_coef_value);
	return py_x_list;
}

/* Module method table */
static PyMethodDef optm_Methods[] = {
	{"fmin_lp", fmin_lp, METH_VARARGS, "Linear Programming"},
	{ NULL, NULL, 0, NULL}
};

/* Module structure */
static struct PyModuleDef optm_module = {
	PyModuleDef_HEAD_INIT,
	"clib_optm",
	"Optimization Module",
	-1,
	optm_Methods
};

/* Module initialization function */
PyMODINIT_FUNC PyInit_clib_optm(void) {
	return PyModule_Create(&optm_module);
}
