"""
Implements a class for performing Basis Adaptation using the Active Subspace method.

Author: Panos Tsilifis
Date: 2/7/2017

"""

__all__ = ['ActiveSubspaceAdaptation', 'GaussianAdaptation', 'QuadraticAdaptation']


import numpy as np
import scipy.stats as st
from . import PolyBasis
from . import MonicPoly
from . import Hermite1d
from . import Legendre1d
from . import Laguerre1d
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

    def __init__(self, num_dim, pol_type = 'H', name = 'Active Subspace Basis Adaptation'):
        """
        Initialize the object
        """
        super(ActiveSubspaceAdaptation, self).__init__(num_dim, pol_type = pol_type, name = name)

    def _stiffness_K(self, deg, i, j):
        """
        Computes the stiffness matrix K^{ij} with entries 
        K_{ab}^{ij} = E[\frac{d psi_a (xi)}{d xi_i} \frac{d psi_b(xi)}{d xi_j}]
        """
        assert isinstance(deg, int)
        assert isinstance(i, int)
        assert isinstance(j, int)
    	assert i < self._inp_dim and j < self._inp_dim
        #rvs = [st.norm()] * self._inp_dim
        
        #pol = orthpol.ProductBasis(rvs, degree = deg)
        pol = PolyBasis(self._inp_dim, deg)
        terms = pol.mi_terms(self._inp_dim, deg)
        Q = terms.shape[0]
        stiff = np.zeros((Q,Q))
        if self._poly_type == 'Legendre' or 'L':
            for k in range(Q):
                for l in range(Q):
                    alpha = terms[k,:]
                    beta = terms[l,:]
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
        elif self._poly_type == 'Hermite' or 'H':
            for k in range(Q):
                for l in range(Q):
                    alpha = terms[k,:]
                    beta = terms[l,:]
                    if alpha[i] - 1 < 0 or beta[j] - 1 < 0:
                        stiff[k,l] = 0
                    else:
                        C = np.sqrt(alpha[i] * beta[j])
                        alpha[i] = alpha[i] - 1
                        beta[j] = beta[j] - 1
                        stiff[k,l] = C * diracND(alpha, beta)
        return stiff

    def _grad_covar(self, deg, coeffs):
        """
        Computes the covariance matric of the gradient vector 
        of the QoI given a set of coefficient values "coeffs"
        corresponding to the currently availble full dimensional
        PCE of the QoI.
        """
        rvs = None
        #if self._poly_type == 'Hermite' or 'H':
        #    rvs = [st.norm()] * self._inp_dim
        pol = PolyBasis(self._inp_dim, deg, self._poly_type[0])
        #elif self._poly_type == 'Legendre' or 'L':
        #    rvs = [st.uniform()] * self._inp_dim
        #pol = orthpol.ProductBasis(rvs, degree = deg)
        terms = pol.mi_terms(self._inp_dim, deg)
        Q = terms.shape[0]
        assert Q == coeffs.shape[0]
        Grad_lo = np.zeros((self._inp_dim, self._inp_dim))
        for i in range(self._inp_dim):
            for j in range(i+1):
                stiff = self._stiffness_K(deg, i, j)
                Grad_lo[i,j] = np.dot(np.dot(coeffs.reshape(1,Q), stiff), coeffs.reshape(Q,1))[0,0]
        return Grad_lo + Grad_lo.T - np.diag(np.diag(Grad_lo))

    def rotation(self, coeffs):
        """
        Computes the rotation matrix.
        """
        raise NotImplementedError('Only my children compute the rotation matrix directly!')


class GaussianAdaptation(ActiveSubspaceAdaptation):
    """
    A class that represents a Polynomial Chaos expansion with adapted basis using the 
    Gaussian adaptation method. 
    """
    def __init__(self, num_dim, name = 'Gaussian Basis Adaptation'):
        super(GaussianAdaptation, self).__init__(num_dim, pol_type = 'H', name = name)

    def rotation(self, coeffs):
        assert coeffs.shape[0] > self._inp_dim
        coeffs = coeffs[:self._inp_dim + 1]
        C = self._grad_covar(1, coeffs)
        [l, v] = np.linalg.eigh(C)
        return v[:,::-1].T


class QuadraticAdaptation(ActiveSubspaceAdaptation):
    """
    A class that represents a Polynomial Chaos expansion with adapted basis using the 
    Quadratic adaptation method. 
    """
    def __init__(self, num_dim, name = 'Quadratic Basis Adaptation'):
        super(QuadraticAdaptation, self).__init__(num_dim, pol_type = 'H', name = name)

    def rotation(self, coeffs):
        Q = self._inp_dim + 1 + self._inp_dim * (self._inp_dim + 1) / 2
        assert coeffs.shape[0] >= Q
        coeffs = coeffs[:Q]
        coeffs[:self._inp_dim + 1] = 0
        C = self._grad_covar(2, coeffs)
        [l, v] = np.linalg.eigh(C)
        return v[:,::-1].T
        
