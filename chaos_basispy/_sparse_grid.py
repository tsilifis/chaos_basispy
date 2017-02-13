"""
Class that generates quadrature points corresponding to sparse grids.
Requirement: Dakota (https://dakota.sandia.gov/)

Author: Panagiotis Tsilifis
Date: 2/10/2017
"""

__all__ = ['SparseGrid']

from os import system
import numpy as np

class SparseGrid(object):
	"""

	"""

	# Name of the SparseGrid object.
	__name__ = None

	# The level of the quadrature rule.
	_level = None 

	# The type of variables.
	_variables = None

	# The parameters of the variables. 
	_params = None

	@property
	def level(self):
		return self._level

	@level.setter
	def level(self, value):
		self._level = value

	@property
	def variables(self):
		return self._type

	@variables.setter
	def variables(self, value):
		self._type = value

	@property
	def params(self):
		return self._params

	@params.setter
	def params(self, value):
		assert value.shape[0] == 2
		if self._variables == "N" or "Normal":
			assert value[1] > 0
			self._params = value
		else:
			self._params = value

	def __init__(self, level, variables, params = None, name = "Sparse Grid Generator"):

		assert isinstance(level, int)
		self._level = level
		if variables == "N" or "Normal":
			self._variables = variables
			if params is not None:
				assert params.shape[0] == 2
				assert params[1] > 0
				self._params = params
			else:
				self._params = np.array([0., 1.])
		elif variables == "U" or "Uniform":
			self._variables = variables
			if params is not None:
				assert params.shape[0] == 2
				self._params = params
			else:
				self._params = np.array([-1., 1.])
		else:
			raise RuntimeError('Currently only Normal and Uniform grids are supported !\n Type should be "Normal" ("N") or "Uniform" ("U").')
		self.__name__ = str(name)

	def generate_grid(self, dim):
		"""
		Generates a sparse grid of dimension "dim" with default parameters. 
		"""
		assert isinstance(dim, int)
		name = 'sparse_grid_gen.in'
		f = open(name, 'w')
		if self._variables == "N" or "Normal":
			text = ['method\n', '	polynomial_chaos\n', '	sparse_grid_level=3\n', '	output verbose\n', 
					'variables\n', '	normal_uncertain ='+str(dim)+'\n', '	means ='+str(dim)+'*'+str(self._params[0])+'\n', '	std_deviations ='+str(dim)+'*'+str(self._params[1])+'\n',
					'interface\n', '	direct\n', "	analysis_driver = 'rosenbrock'\n", 'responses\n', 
					"	response_functions = 1\n", "	no_gradients\n", "	no_hessians\n"]
		elif self._variables == "U" or "Uniform":
			text = ['method\n', '	polynomial_chaos\n', '	sparse_grid_level=3\n', '	output verbose\n', 
					'variables\n', '	uniform_uncertain ='+str(dim)+'\n', '	lower_bounds ='+str(dim)+'* '+str(self._params[0])+'\n', '	upper_bounds ='+str(dim)+'*'+str(self._params[1])+'\n',
					'interface\n', '	direct\n', "	analysis_driver = 'rosenbrock'\n", 'responses\n', 
					"	response_functions = 1\n", "	no_gradients\n", "	no_hessians\n"]

		f.writelines(text)
		f.close()
		system('dakota -i ' + str(f.name))

		out = open('dakota_sparse_tabular.dat', 'rU')
		N = -1
		for line in f:
			N = N + 1
		weights = np.zeros(N)
		grid = np.zeros((N, dim))
		f.readline()

		pos = []
		for i in range(dim):
			pos += [25+i*25]

		for i in range(N):
			line = f.readline()
			print line[:6]
			weights[i] = float(line[6:25])
			grid[i,:] = np.array([float(line[j:j+15]) for j in pos])

		return grid, weights
			

	def __str__(self):
		"""
		Return a string representation of the object.
		"""
		s = 'Name: ' + self.__name__ + '\n'
		s += 'Level: ' + str(self._level) + '\n'
		s += 'Variables: ' + str(self._variables) + '\n'
		return s


