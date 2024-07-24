from typing import Literal

from . import clib_optm

class LinearObjective:
	coef: list[float]
	type: str

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

class LinearConstraint:
	coef: list[float]
	type: str
	type_code: int
	rhs: float

	def __init__(self, coef: list[float], type: Literal['>=', '<=', '=='], rhs: float) -> None:
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
			raise ValueError('Expected either ">=", "<=" or "==", get "{}".'.format(type))
		self.rhs = float(rhs)
	def __repr__(self) -> str:
		return str(self.coef) + ' ' + self.type + ' ' + str(self.rhs)

class LinearProgramming:
	obj: LinearObjective
	consts: list[LinearConstraint]

	def __init__(self, obj: LinearObjective, consts: list[LinearConstraint]) -> None:
		self.obj = obj
		self.consts = consts
	def __repr__(self) -> str:
		re = self.obj.__repr__() + '\n' + 's.t.' + '\n'
		for cons in self.consts:
			re += '     ' + cons.__repr__() + '\n'
		return re
	def solve(self, method: Literal['bland', ] = 'bland', max_iter: int = 1000):
		# obj_coef
		if self.obj.type == 'min':
			obj_coef = self.obj.coef
		else:
			obj_coef = [-e for e in self.obj.coef]
		consts_coef = []
		consts_rhs = []
		consts_type = []
		for constraint in self.consts:
			consts_coef.append(constraint.coef)
			consts_rhs.append(constraint.rhs)
			consts_type.append(constraint.type_code)
		m = len(self.consts)
		n = len(obj_coef)
		return clib_optm.fmin_lp(
			m, n, max_iter, method,
			obj_coef, consts_coef, consts_rhs, consts_type
		)
