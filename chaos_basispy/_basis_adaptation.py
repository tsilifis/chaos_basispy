"""
Implements a class representing an adapted Polynomial Chaos expansion.

Author : Panos Tsilifis
Date : 2/7/2017

"""

__all__ = ['BasisAdaptation']


import numpy as np 


class BasisAdaptation(object):
	"""
	A class that represents a Polynomial Chaos expansion 
	expressed in terms of a rotated basis. 

	"""

	# A name for the expansion
	__name__ = None

	# The type of polynomials in the expansion
	_poly_type = None

	# The initial input dimension.
	_inp_dim = None

	# The initial highest order of the Chaos space.
	_chaos_order = None

	# The series coefficients with respect to the initial input xi.
	_chaos_coeffs = None

	# The number of coefficients with respect to the initial input xi.
	_num_chaos_coeffs = None

	@property
	def inp_dim(self):
		"""
		:getter: The initial input dimension
		"""
		return self._inp_dim

	@inp_dim.setter
	def inp_dim(self, value):
		"""
		Set the input dimension.
		"""
		self._inp_dim = value

	@property
	def chaos_order(self):
		"""
		:getter: The initial highest order of the Chaos space.
		"""
		return self._chaos_order

	@chaos_order.setter
	def chaos_order(self, value):
		"""
		Set the initial highest order of the Chaos space.
		"""
		self._chaos_order = value

	@property
	def num_chaos_coeffs(self):
		"""
		:getter: The number of coefficients with respect to the initial input xi.
		"""
		return self._num_chaos_coeffs

	@num_chaos_coeffs.setter
	def num_chaos_coeffs(self, value):
		"""
		Set the number of coefficients with respect to the initial input xi.
		"""
		self._num_chaos_coeffs = value

	@property
	def poly_type(self):
		"""
		:getter: The type of polynomials in the expansion 
		(currently supports only Hermite and Legendre).
		"""
		return self._poly_type

	@poly_type.setter
	def poly_type(self, value):
		"""
		Set the type of polynomials to be either Hermite of Legendre.
		"""
		if value == 'Hermite' or "H" or 'Legendre' or "L":
			self._poly_type = value
		else:
			raise RuntimeError('The polynomials should be either Hermite of Legendre!')

	def __init__(self, num_dim, num_chaos_coeffs = None, chaos_order = None, chaos_coeffs = None, pol_type = 'Hermite', name = 'Adapted PC expansion'):
		"""
		Initialize the object.
		"""
		assert isinstance(num_dim, int)
		self._inp_dim = num_dim
		if num_chaos_coeffs is not None:
			self._num_chaos_coeffs = num_chaos_coeffs
		if chaos_order is not None:
			self._chaos_order = chaos_order
		if chaos_coeffs is not None:
			self._chaos_coeffs = chaos_coeffs
		if pol_type == 'Hermite' or 'Legendre' or 'H' or 'L':
			self._poly_type = pol_type
		else:
			raise RuntimeError('The polynomials should be either Hermite of Legendre!')
		self.__name__ = str(name)


	def __str__(self):
		"""
		Return a string representation of the object.
		"""
		s = 'Name: ' + self.__name__ + '\n'
		s += 'Type of polynomials: ' + self._poly_type + '\n'
		s += 'Initial input dimension: ' + str(self._inp_dim) + '\n'
		s += 'Order of expansion: ' + str(self._chaos_order) + '\n'
		return s








