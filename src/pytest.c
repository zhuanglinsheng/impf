#define PY_SSIZE_T_CLEAN  /* Make "s#" use Py_ssize_t rather than int. */
#include <Python.h>

static PyObject* test(PyObject* self, PyObject* args) {
	return Py_BuildValue("s", "Greeting from test module!");
}

/* Module method table */
static PyMethodDef test_Methods[] = {
	{"test", test, METH_VARARGS, "greeting"},
	{ NULL, NULL, 0, NULL}
};

/* Module structure */
static struct PyModuleDef test_module = {
	PyModuleDef_HEAD_INIT,
	"clib_test",
	"Test Module",
	-1,
	test_Methods
};

/* Module initialization function */
PyMODINIT_FUNC PyInit_clib_test(void) {
	return PyModule_Create(&test_module);
}
