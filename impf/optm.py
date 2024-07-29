from typing import Literal
from . import _clib_optm
from . import _io

class LinearObjective:
	coef: list[float]
	type: Literal['max', 'min']

	def __init__(self, coef: list[float], type: Literal['max', 'min']) -> None:
		self.coef = [float(e) for e in coef]
		if type == 'max':
			self.type = 'max'
		elif type == 'min':
			self.type = 'min'
		else:
			raise ValueError('Expected either "max" or "min", get "{}".'.format(type))
	def __repr__(self) -> str:
		return self.type + '  ' + str(self.coef)
	def get_value(self, x: list[float]) -> float:
		val = 0
		for (x_i, c_i) in zip(x, self.coef):
			val += x_i * c_i
		return val


class LinearConstraint:
	coef: list[float]
	type: Literal['==', '>=', '<=']
	type_code: Literal[0, 1, 2]
	rhs: float

	def __init__(self, coef: list[float], type: Literal['==', '>=', '<='], rhs: float) -> None:
		self.coef = [float(e) for e in coef]
		if type == '==':
			self.type = '=='
			self.type_code = 0
		elif type == '>=':
			self.type = '>='
			self.type_code = 1
		elif type == '<=':
			self.type = '<='
			self.type_code = 2
		else:
			raise ValueError('Expected either "==", ">=" or "<=", get "{}".'.format(type))
		self.rhs = float(rhs)
	def __repr__(self) -> str:
		return str(self.coef) + ' ' + self.type + ' ' + str(self.rhs)

class LinearProgramming:
	obj: LinearObjective
	consts: list[LinearConstraint]
	bounds: list[tuple[float, float]]

	def __init__(self, obj: LinearObjective,
			consts: list[LinearConstraint],
			bounds: list[tuple[float, float]] = None) -> None:
		if (bounds is not None) and len(bounds) != len(obj.coef):
			raise ValueError('len(obj.coef) = {} but len(bounds) = {}'.format(len(obj.coef), len(bounds)))
		self.obj = obj
		self.consts = consts
		self.bounds = bounds
	def __repr__(self) -> str:
		blk = '     '
		re = self.obj.__repr__() + '\n' + 's.t.' + '\n'
		for cons in self.consts:
			re += blk + cons.__repr__() + '\n'
		if self.bounds is None:
			re += blk + 'all variables >= 0\n'
		else:
			for i in range(len(self.bounds)):
				re += blk + str(self.bounds[i][0]) + ' <= x' + str(i) + \
					' <= ' + str(self.bounds[i][1]) + '\n'
		return re

	def readMPS(file: str):
		tokens = _io.readMPS(file)
		obj = LinearObjective(tokens['obj_coef'], 'min')
		constraints = []
		con_coef = tokens['constraints_coef']
		con_type = tokens['constraints_type']
		con_rhs = tokens['constraints_rhs']
		for (coef, type, rhs) in zip(con_coef, con_type, con_rhs):
			constraints.append(LinearConstraint(coef, type, rhs))
		return LinearProgramming(obj, constraints, tokens['vars_bound'])

	def solve(self, method: Literal['dantzig'] = 'dantzig', max_iter: int = 1000) -> dict:
		if self.obj.type == 'min':
			obj_coef = self.obj.coef
		else:
			obj_coef = [-e for e in self.obj.coef]
		m = len(self.consts)
		n = len(obj_coef)
		consts_coef = []
		consts_rhs = []
		consts_type = []
		for constraint in self.consts:
			consts_coef.append(constraint.coef)
			consts_rhs.append(constraint.rhs)
			consts_type.append(constraint.type_code)
		re = dict()
		x, code = _clib_optm.wrapper_impf_lp_simplex(
			m, n, max_iter, method,
			self.bounds, obj_coef, consts_coef, consts_rhs, consts_type
		)
		re['x'] = x

		if code == 0:
			re['state'] = 'Success'
		elif code == 1:
			re['state'] = 'MemoryAllocError'
		elif code == 2:
			re['state'] = 'CondUnsatisfied'
		elif code == 3:
			re['state'] = 'ExceedIterLimit'
		elif code == 4:
			re['state'] = 'Singularity'
		elif code == 5:
			re['state'] = 'OverDetermination'
		elif code == 6:
			re['state'] = 'Unboundedness'
		elif code == 7:
			re['state'] = 'Infeasibility'
		elif code == 8:
			re['state'] = 'Degeneracy'
		elif code == 9:
			re['state'] = 'PrecisionError'
		return re
