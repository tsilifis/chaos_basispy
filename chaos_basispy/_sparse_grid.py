"""
Class that generates quadrature points corresponding to sparse grids.
Requirement: Dakota (https://dakota.sandia.gov/)

Author: Panagiotis Tsilifis
Date: 2/10/2017
"""

__all__ = ['SparseGrid']

from os import system

class SparseGrid(object):
	"""

	"""
	__name__ = None

	_level = None 

	_variables = None

	@property
	def level(self):
		return self._level

	@level.setter
	def level(self, value):
		self._level = value

	@property
	def variables(self):
		return self._type

	@type.setter
	def variables(self, value):
		self._type = value

	def __init__(self, level, variables, name = "Sparse Grid Generator"):

		assert isinstance(leve, int)
		self._level = level
		if type == "N" or "Normal":
			self._variables = variables
		elif type == "U" or "Uniform":
			se;f._variables = variables
		else:
			raise RuntimeError('Currently only Normal and Uniform grids are supported !\n Type should be "Normal" ("N") or "Uniform" ("U").')
		self.__name__ = str(name)

	def generate_grid(dim):
		"""
		Generates a sparse grid of dimension "dim" with default parameters. 
		"""
		assert isinstance(dim, int)
		name = 'sparse_grid_gen.in'
		f = open(name, 'w')
		if self._variables == "N" or "Normal":
			text = ['method\n', '	polynomial_chaos\n', '	sparse_grid_level=3\n', '	output verbose\n', 
					'variables\n', '	normal_uncertain ='+str(dim)+'\n', '	means ='+str(dim)+'*0\n', '	std_deviations ='+str(dim)+'*1\n',
					'interface\n', '	direct\n', "	analysis_driver = 'rosenbrock'\n", 'responses\n', 
					"	response_functions = 1\n", "	no_gradients\n", "	no_hessians\n"]
		elif self._variables == "U" or "Uniform":
			text = ['method\n', '	polynomial_chaos\n', '	sparse_grid_level=3\n', '	output verbose\n', 
					'variables\n', '	uniform_uncertain ='+str(dim)+'\n', '	lower_bounds ='+str(dim)+'* -1.\n', '	upper_bounds ='+str(dim)+'*1\n',
					'interface\n', '	direct\n', "	analysis_driver = 'rosenbrock'\n", 'responses\n', 
					"	response_functions = 1\n", "	no_gradients\n", "	no_hessians\n"]

		f.writelines(text)
		f.close()
		system('dakota -i ' + str(f.name))

	def __str__(self):
		"""
		Return a string representation of the object.
		"""
		s = 'Name: ' + self.__name__ + '\n'
		s += 'Level: ' + self._level + '\n'
		s += 'Variables: ' + str(self._variables) + '\n'
		return s


