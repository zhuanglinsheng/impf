from . import expr
from . import optm
from typing import Literal
import math

MPS_Mode_t = Literal['', 'ROWS', 'COLUMNS', 'RHS', 'RANGES', 'BOUNDS']
Bound_t = Literal['LO', 'UP', 'FX', 'FR', 'MI', 'PL', 'BV', 'LI', 'UI', 'SC']
Constraint_t = Literal['==', '>=', '<=']

def readMPS(filepath: str) -> dict:
	'''
	Read tokens of MPS file.
	Return
		dict {
			'prob_name'       : str,
			'vars_pool'       : list[str],
			'vars_bound'      : Option[list[tuple[float, float]]],
			'vars_bound_type' : Option[list[Bound_t]],
			'obj_name'        : list[str],
			'obj_coef'        : list[float],
			'constraints_name': list[str],
			'constraints_type': list[Constraint_t],
			'constraints_coef': list[list[float]],
			'constraints_rhs' : list[float]
		}
	where
		Bound_t = Literal['LO', 'UP', 'FX', 'FR', 'MI', 'PL', 'BV', 'LI', 'UI', 'SC']
		Constraint_t = Literal['==', '>=', '<=']
	Reference: https://lpsolve.sourceforge.net/5.5/mps-format.htm
	'''
	re = dict()
	re['prob_name'] = ''
	re['vars_pool'] = []
	re['vars_bound'] = None
	re['vars_bound_type'] = None
	re['obj_name'] = ''
	re['obj_coef'] = []
	re['constraints_name'] = []
	re['constraints_type'] = []
	re['constraints_coef'] = []
	re['constraints_rhs'] = []

	mode: MPS_Mode_t = ''
	num_vars = 0
	num_cons = 0
	mps_has_bound = False

	# first loop: get information
	with open(filepath, 'r') as file:
		for line in file:
			line_trimmed = line.strip()

			if line_trimmed[:4] == 'NAME':
				mode = ''
				continue
			elif line_trimmed == 'ROWS':
				mode = 'ROWS'
				continue
			elif line_trimmed == 'COLUMNS':
				mode = 'COLUMNS'
				continue
			elif line_trimmed == 'RHS':
				mode = 'RHS'
				continue
			elif line_trimmed == 'RANGES':
				mode = 'RANGES'
				continue
			elif line_trimmed == 'BOUNDS':
				mode = 'BOUNDS'
				mps_has_bound = True
				continue
			elif line_trimmed == 'ENDATA':
				break

			if mode == 'ROWS':
				row_type, row_name = line_trimmed.split()

				if row_type == 'N':
					re['obj_name'] = row_name
				elif row_type == 'E':
					re['constraints_type'].append('==')
					re['constraints_name'].append(row_name)
				elif row_type == 'G':
					re['constraints_type'].append('>=')
					re['constraints_name'].append(row_name)
				elif row_type == 'L':
					re['constraints_type'].append('<=')
					re['constraints_name'].append(row_name)
			if mode == 'COLUMNS':
				var_infos = line_trimmed.split()
				var_name = var_infos[0]

				if var_name not in re['vars_pool']:
					re['vars_pool'].append(var_name)
		num_vars = len(re['vars_pool'])
		num_cons = len(re['constraints_name'])
		zeros = [0.] * num_vars
		re['obj_coef'] = zeros.copy()
		re['constraints_coef'] = [zeros.copy() for _ in range(num_cons)]
		re['constraints_rhs'] = [0.] * num_cons
		if mps_has_bound:
			re['vars_bound'] = [(0., math.inf) for _ in range(num_vars)]

	# second loop: fill in value
	with open(filepath, 'r') as file:
		for line in file:
			line_trimmed = line.strip()

			if line_trimmed[:4] == 'NAME':
				mode = ''
				re['prob_name'] = line_trimmed[4:].strip()
				continue
			elif line_trimmed == 'ROWS':
				mode = 'ROWS'
				continue
			elif line_trimmed == 'COLUMNS':
				mode = 'COLUMNS'
				continue
			elif line_trimmed == 'RHS':
				mode = 'RHS'
				continue
			elif line_trimmed == 'RANGES':
				mode = 'RANGES'
				continue
			elif line_trimmed == 'BOUNDS':
				mode = 'BOUNDS'
				continue
			elif line_trimmed == 'ENDATA':
				break

			if mode == 'COLUMNS':
				var_infos = line_trimmed.split()
				var_name = var_infos[0]
				var_idx = re['vars_pool'].index(var_name)

				for i in range(1, len(var_infos) - 1, 2):
					var_belong = var_infos[i]
					var_coef = float(var_infos[i + 1])
					if var_belong == re['obj_name']:
						re['obj_coef'][var_idx] = var_coef
					else:
						var_belong_idx = re['constraints_name'].index(var_belong)
						re['constraints_coef'][var_belong_idx][var_idx] = var_coef
			if mode == 'RHS':
				rhs_infos = line_trimmed.split()
				for i in range(1, len(rhs_infos) - 1, 2):
					rhs_belong = rhs_infos[i]
					rhs_value = float(rhs_infos[i + 1])
					rhs_belong_idx = re['constraints_name'].index(rhs_belong)
					re['constraints_rhs'][rhs_belong_idx] = rhs_value
			if mode == 'RANGES':
				pass
			if mode == 'BOUNDS':
				pass
	return re

def parseMPS(tokens: dict) -> optm.LinearProgramming:

	obj = optm.LinearObjective(tokens['obj_coef'], 'min')
	constraints = []
	con_coef = tokens['constraints_coef']
	con_type = tokens['constraints_type']
	con_rhs = tokens['constraints_rhs']

	for (coef, type, rhs) in zip(con_coef, con_type, con_rhs):
		constraints.append(optm.LinearConstraint(coef, type, rhs))

	return optm.LinearProgramming(obj, constraints, tokens['vars_bound'])

