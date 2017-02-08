"""
Implements a class for performing Basis Adaptation using the Active Subspace method.

Author: Panos Tsilifis
Date: 2/7/2017

"""

__all__ = ['AS_adaptation']


import numpy as np
from . import BasisAdaptation

class ActiveSubspaceAdaptation(BasisAdaptation):
	"""
	A class the represents a Polynomial Chaos expansion with adapted basis using the 
	Active Subspace method. 

	"""

	def __init__(self, num_dim, name = 'Active Subspace Basis Adaptation'):
		"""
		Initialize the object
		"""
		super(ActiveSubspaceAdaptation, self).__init__(num_dim, name = name)

	