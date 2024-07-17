import ctypes
import ctypes.util

# load impf library
lib_impf_path = ctypes.util.find_library('impf')
if lib_impf_path is None:
	raise OSError('Error loading dynamic library: impf')
_lib_impf: ctypes.CDLL = ctypes.CDLL(lib_impf_path)

def get_lib_impf_instance() -> ctypes.CDLL:
	return _lib_impf
