"""
Implements a class for performing Basis Adaptation using the Active Subspace method.

Author: Panos Tsilifis
Date: 2/7/2017

"""

__all__ = ['ActiveSubspaceAdaptation', 'GaussianAdaptation', 'QuadraticAdaptation']


import numpy as np
import scipy.stats as st
import orthpol
from . import BasisAdaptation
from . import dirac1D
from . import diracND
from . import grad_inner_prod_Legendre
from . import gradgrad_inner_prod_Legendre

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

    def _stiffness_K(self, deg, i, j):
        """
        Computes the stiffness matrix K^{ij} with entries 
        K_{ab}^{ij} = E[\frac{d psi_a (xi)}{d xi_i} \frac{d psi_b(xi)}{d xi_j}]
        """
        assert isinstance(deg, int)
        assert isinstance(i, int)
        assert isinstance(j, int)
    	assert i < self._inp_dim and j < self._inp_dim
        rvs = [st.norm()] * self._inp_dim
        pol = orthpol.ProductBasis(rvs, degree = deg)
        Q = len(pol._terms)
        stiff = np.zeros((Q,Q))
        if self._poly_type == 'Legendre':
            for k in range(Q):
                for l in range(Q):
                    alpha = pol._terms[k]
                    beta = pol._terms[l]
                    if i == j :
                        a_m = np.delete(alpha, i)
                        b_m = np.delete(beta, i)
                        a_i = alpha[i]
                        b_i = beta[i]
                        E_ii = gradgrad_inner_prod_Legendre(a_i,b_i)
                        stiff[k,l] = diracND(a_m, b_m) * E_ii
                    else:
                        a_m = np.delete(alpha, [i,j])
                        b_m = np.delete(beta, [i,j])
                        a_i = alpha[i]
                        a_j = alpha[j]
                        b_i = beta[i]
                        b_j = beta[j]
                        E_i = grad_inner_prod_Legendre(a_i, b_i)
                        E_j = grad_inner_prod_Legendre(b_j, a_j)
                        stiff[k,l] = diracND(a_m, b_m) * E_i * E_j
        elif self._poly_type == 'Hermite':
            for k in range(Q):
                for l in range(Q):
                    alpha = pol._terms[k]
                    beta = pol._terms[l]
                    if alpha[i] - 1 < 0 or beta[j] - 1 < 0:
                        stiff[k,l] = 0
                    else:
                        C = np.sqrt(alpha[i] * beta[j])
                        alpha[i] = alpha[i] - 1
                        beta[j] = beta[j] - 1
                        stiff[k,l] = C * diracND(alpha, beta)
        return stiff


class GaussianAdaptation(ActiveSubspaceAdaptation):
    """
    docstring for GaussianAdaptation"ActiveSubspaceAdaptation 
    """
    def __init__(self, num_dim, name = 'Gaussian Basis Adaptation'):
        super(GaussianAdaptation, self).__init__(num_dim, name = name)

class QuadraticAdaptation(ActiveSubspaceAdaptation):
    """
    """
    def __init__(self, num_dim, name = 'Quadratic Basis Adaptation'):
        super(QuadraticAdaptation, self).__init__(num_dim, name = name)
        
