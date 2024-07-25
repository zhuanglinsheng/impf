from __future__ import annotations
from typing import Literal, Union
import math


class LinearExpression:
	left: LinearExpression
	right: LinearExpression
	type: Literal['add', 'sub', 'mul']
	degree: int

	def __init__(self, left: LinearExpression, right: LinearExpression,
			type: Literal['add', 'sub', 'mul']) -> None:
		self.left = left
		self.right = right
		self.type = type
		self.degree = 0
		if type == 'add' or type == 'sub':
			self.degree = max(left.degree, right.degree)
		elif type == 'mul':
			self.degree = left.degree + right.degree
		if self.degree > 1:
			raise ValueError('Not a Linear Expression')
	def __repr__(self) -> str:
		if self.type == 'add':
			type_str = '+'
		elif self.type == 'sub':
			type_str = '-'
		elif self.type == 'mul':
			type_str = '*'
		re = ''
		if self.left._is_atom():
			re += str(self.left)
		else:
			re += '(' + str(self.left) + ')'
		re += ' ' + type_str + ' '
		if self.right._is_atom():
			re += str(self.right)
		else:
			re += '(' + str(self.right) + ')'
		return re

	def tranform(self, var_pool: list[ScalarVariable]) -> tuple[list[float], float]:
		'''return (coefficients, constant)
		'''
		coefficients = [0.0] * len(var_pool)
		constant = 0
		if self.left._is_atom() and (self.type == 'add' or self.type == 'sub'):
			if isinstance(self.left, ScalarConstant):
				constant += self.left.data
			elif isinstance(self.left, ScalarVariable):
				idx = var_pool.index(self.left)
				coefficients[idx] = 1. if self.type == 'add' else -1.0

	def _is_atom(self) -> bool:
		return self.type is None and self.left is None and self.right is None

	def __add__(self, other: Union[int, float, LinearExpression]) -> LinearExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return LinearExpression(self, other_expr, 'add')
		elif isinstance(LinearExpression):
			return LinearExpression(self, other, 'add')
		else:
			raise ValueError('Expect Union[int, float, LinearExpression], get {}'.format(type(other)))
	def __radd__(self, other: Union[int, float, LinearExpression]) -> LinearExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return LinearExpression(other_expr, self, 'add')
		elif isinstance(LinearExpression):
			return LinearExpression(other, self, 'add')
		else:
			raise ValueError('Expect Union[int, float, LinearExpression], get {}'.format(type(other)))
	def __sub__(self, other: Union[int, float, LinearExpression]) -> LinearExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return LinearExpression(self, other_expr, 'sub')
		elif isinstance(LinearExpression):
			return LinearExpression(self, other, 'sub')
		else:
			raise ValueError('Expect Union[int, float, LinearExpression], get {}'.format(type(other)))
	def __mul__(self, other: Union[int, float, LinearExpression]) -> LinearExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return LinearExpression(self, other_expr, 'mul')
		elif isinstance(LinearExpression):
			return LinearExpression(self, other, 'mul')
		else:
			raise ValueError('Expect Union[int, float, LinearExpression], get {}'.format(type(other)))



class ScalarConstant(LinearExpression):
	data: float
	degree: int

	def __init__(self, data: float) -> None:
		super().__init__(None, None, None)
		self.data = data
		self.degree = 0
	def __repr__(self) -> str:
		return str(self.data)
	def _is_atom(self) -> bool:
		return True


class ScalarVariable(LinearExpression):
	type: Literal['real', 'integer', 'binary']
	name: str
	data: float
	lb: float
	ub: float
	degree: int

	def __init__(self, name: str = '', type: Literal['real', 'integer', 'binary'] = 'real',
			lb: float = -math.inf, ub: float = math.inf, data: float = 0.0) -> None:
		super().__init__(None, None, None)
		self.name = name
		self.type = type
		self.data = data
		self.lb = lb
		self.ub = ub
		self.degree = 1
	def __repr__(self) -> str:
		if len(self.name) > 0:
			return self.name
		else:
			return 'var'
	def _is_atom(self) -> bool:
		return True
