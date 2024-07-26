from __future__ import annotations
from typing import Literal, Union
import math


class Expression:
	def __init__(self) -> None:
		pass
	def _is_atom(self) -> bool:
		return False


class BinExpression(Expression):
	left: BinExpression
	right: BinExpression
	type: Literal['add', 'sub', 'mul']

	def __init__(self, left: BinExpression, right: BinExpression,
			type: Literal['add', 'sub', 'mul']) -> None:
		self.left = left
		self.right = right
		self.type = type
	def __repr__(self) -> str:
		if self.type == 'add':
			type_str = '+'
		elif self.type == 'sub':
			type_str = '-'
		elif self.type == 'mul':
			type_str = '*'
		re = ''
		if self.type == 'mul' and not self.left._is_atom():
			re += '(' + str(self.left) + ')'
		else:
			re += str(self.left)
		re += ' ' + type_str + ' '
		if self.type == 'mul' and not self.right._is_atom():
			re += '(' + str(self.right) + ')'
		else:
			re += str(self.right)
		return re
	def get_degree(self) -> int:
		if self.type == 'add' or self.type == 'sub':
			return max(self.left.get_degree(), self.right.get_degree())
		if self.type == 'mul':
			return self.left.get_degree() + self.right.get_degree()
	def __add__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(self, other_expr, 'add')
		elif isinstance(other, BinExpression):
			return BinExpression(self, other, 'add')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))
	def __radd__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(other_expr, self, 'add')
		elif isinstance(other, BinExpression):
			return BinExpression(other, self, 'add')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))
	def __sub__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(self, other_expr, 'sub')
		elif isinstance(other, BinExpression):
			return BinExpression(self, other, 'sub')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))
	def __rsub__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(other_expr, self, 'sub')
		elif isinstance(other, BinExpression):
			return BinExpression(other, self, 'sub')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))
	def __mul__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(self, other_expr, 'mul')
		elif isinstance(other, BinExpression):
			return BinExpression(self, other, 'mul')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))
	def __rmul__(self, other: Union[int, float, BinExpression]) -> BinExpression:
		if isinstance(other, int) or isinstance(other, float):
			other_expr = ScalarConstant(float(other))
			return BinExpression(other_expr, self, 'mul')
		elif isinstance(other, BinExpression):
			return BinExpression(other, self, 'mul')
		else:
			raise ValueError('Expect Union[int, float, BinExpression], get {}'.format(type(other)))


class LinearExpression(BinExpression):
	def __init__(self, expr: BinExpression) -> None:
		super().__init__(expr.left, expr.right, expr.type)
		if self.get_degree() > 1:
			raise ValueError('"expr" is not a linear expression: {}'.format(str(expr)))

	def tranform(self, var_pool: list[ScalarVariable]) -> tuple[list[float], float]:
		'''return (coefficients: list[float], constant: float)
		'''
		if self.type == 'add':
			if (not self.left._is_atom()) and (not self.right._is_atom()):
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				coefficients = [a + b for (a, b) in zip(left_coef, right_coef)]
				constant = left_const + right_const
			if (not self.left._is_atom()) and self.right._is_atom():
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				if isinstance(self.right, ScalarConstant):
					coefficients = left_coef
					constant = left_const + self.right.data
				elif isinstance(self.right, ScalarVariable):
					coefficients = left_coef
					coefficients[var_pool.index(self.right)] += 1
					constant = left_const
			if self.left._is_atom() and (not self.right._is_atom()):
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				if isinstance(self.left, ScalarConstant):
					coefficients = right_coef
					constant = left_const + self.right.data
				elif isinstance(self.right, ScalarVariable):
					coefficients = right_coef
					coefficients[var_pool.index(self.left)] += 1
					constant = left_const
			if self.left._is_atom() and self.right._is_atom():
				if isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					constant = self.left.data + self.right.data
				elif isinstance(self.left, ScalarVariable) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.left)] = 1.
					constant = self.right.data
				elif isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarVariable):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.right)] = 1.
					constant = self.left.data
		elif self.type == 'sub':
			if (not self.left._is_atom()) and (not self.right._is_atom()):
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				coefficients = [a - b for (a, b) in zip(left_coef, right_coef)]
				constant = left_const - right_const
			if (not self.left._is_atom()) and self.right._is_atom():
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				if isinstance(self.right, ScalarConstant):
					coefficients = left_coef
					constant = left_const - self.right.data
				elif isinstance(self.right, ScalarVariable):
					coefficients = left_coef
					coefficients[var_pool.index(self.right)] -= 1
					constant = left_const
			if self.left._is_atom() and (not self.right._is_atom()):
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				if isinstance(self.left, ScalarConstant):
					coefficients = [-e for e in right_coef]
					constant = left_const - self.right.data
				elif isinstance(self.right, ScalarVariable):
					coefficients = [-e for e in right_coef]
					coefficients[var_pool.index(self.left)] += 1
					constant = -left_const
			if self.left._is_atom() and self.right._is_atom():
				if isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					constant = self.left.data - self.right.data
				elif isinstance(self.left, ScalarVariable) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.left)] = 1.
					constant = -self.right.data
				elif isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarVariable):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.right)] = -1.
					constant = self.left.data
		elif self.type == 'mul':
			if (not self.left._is_atom()) and (not self.right._is_atom()):
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				left_degen = True
				for e in left_coef:
					if e != 0.:
						left_coef = False
				if left_degen:
					coefficients = [left_const * e for e in right_coef]
					constant = left_const * right_const
				else:
					coefficients = [e * right_const for e in left_coef]
					constant = left_const * right_const
			if (not self.left._is_atom()) and self.right._is_atom():
				left_coef, left_const = LinearExpression(self.left).tranform(var_pool)
				if isinstance(self.right, ScalarConstant):
					coefficients = [e * self.right.data for e in left_coef]
					constant = left_const * self.right.data
				elif isinstance(self.right, ScalarVariable):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.right)] = left_const
					constant = 0.
			if self.left._is_atom() and (not self.right._is_atom()):
				right_coef, right_const = LinearExpression(self.right).tranform(var_pool)
				if isinstance(self.left, ScalarConstant):
					coefficients = [self.left.data * e for e in right_coef]
					constant = self.left.data * right_const
				elif isinstance(self.left, ScalarVariable):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.left)] = right_const
					constant = 0.
			if self.left._is_atom() and self.right._is_atom():
				if isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					constant = self.left.data * self.right.data
				elif isinstance(self.left, ScalarVariable) and isinstance(self.right, ScalarConstant):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.left)] = self.right.data
					constant = 0.
				elif isinstance(self.left, ScalarConstant) and isinstance(self.right, ScalarVariable):
					coefficients = [0.] * len(var_pool)
					coefficients[var_pool.index(self.right)] = self.left.data
					constant = 0.
		return coefficients, constant


class ScalarConstant(BinExpression):
	data: float

	def __init__(self, data: float) -> None:
		super().__init__(None, None, None)
		self.data = data
	def __repr__(self) -> str:
		return str(self.data)
	def _is_atom(self) -> bool:
		return True
	def get_degree(self) -> int:
		return 0.


class ScalarVariable(BinExpression):
	type: Literal['real', 'integer', 'binary']
	name: str
	data: float
	lb: float
	ub: float

	def __init__(self, name: str = '', type: Literal['real', 'integer', 'binary'] = 'real',
			lb: float = -math.inf, ub: float = math.inf, data: float = 0.0) -> None:
		super().__init__(None, None, None)
		self.name = name
		self.type = type
		self.data = data
		self.lb = lb
		self.ub = ub
	def __repr__(self) -> str:
		if len(self.name) > 0:
			return self.name
		else:
			return 'var'
	def _is_atom(self) -> bool:
		return True
	def get_degree(self) -> int:
		return 1

