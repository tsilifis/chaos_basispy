""" Basis terms Module 

Author: Panos Tsilifis
Date: 7/10/2017

"""

__all__ = ['PolyBasis', 'MonicPoly', 'Hermite1d', 'Legendre1d', 'Laguerre1d']

import numpy as np
import math
from scipy import misc
import itertools as itls


class Hermite1d(object):
    """
    Class that constructs 1-dimensional Hermite polynomials.
    """


    _degree = None


    def __init__(self, degree = 1):
        """
        Ininializes the object 
        """
        assert isinstance(degree, int)
        assert degree > 0
        self._degree = degree 

    def eval(self, x):
        H = np.zeros(self._degree + 1)
        H[0] = 1.
        H[1] = x
        for i in range(2,H.shape[0]):
            H[i] = (x * H[i-1] - (i-1) * H[i-2] )
        H = H / [math.sqrt(math.factorial(i)) for i in range(H.shape[0])]
        return H

    def __call__(self, x):
        N = x.shape[0]
        H = np.zeros((N, self._degree + 1))
        for i in range(N):
            H[i,:] = self.eval(x[i])
        return H


class Legendre1d(object):
    """
    Class that contains 1-dimensional Legendre polynomials.
    """

    _degree = None

    def __init__(self, degree = 1):
        """
        Initializes the object
        """
        assert isinstance(degree, int)
        assert degree > 0
        self._degree = degree

    def eval(self, x):
        H = np.zeros(self._degree + 1)
        H[0] = 1.
        H[1] = x
        for i in range(2, H.shape[0]):
            H[i] = ( (2.*i-1.)* x * H[i-1] - (i-1.)* H[i-2] ) / i
        H = H / [math.sqrt(1. / (2.*i+1.)) for i in range(H.shape[0])]
        return H

    def __call__(self, x):
        N = x.shape[0]
        H = np.zeros((N, self._degree + 1))
        for i in range(N):
            H[i,:] = self.eval(x[i])
        return H


class Laguerre1d(object):
    """
    Class that contains 1-dimensional Laguerre polynomials.
    """

    _degree = None

    def __init__(self, degree = 1):
        """
        Initializes the object
        """
        assert isinstance(degree, int)
        assert degree > 0
        self._degree = degree

    def eval(self, x):
        H = np.zeros(self._degree + 1)
        H[0] = 1.
        H[1] = 1. - x
        for i in range(2, H.shape[0]):
            H[i] = ((2*(i-1)+1-x)*H[i-1] - (i-1)*H[i-2]) / (i)
        return H

    def __call__(self, x):
        N = x.shape[0]
        H = np.zeros((N, self._degree + 1))
        for i in range(N):
            H[i,:] = self.eval(x[i])
        return H


class MonicPoly(object):
    """
    Class that contains 1-dimensioanl monic polynomials, 
    orthogonal with respect to an arbitrary measure. 
    """

    _degree = None
    _measure = None
    _support = None

    def __init__(self, degree, pdf, supp):
        """
        Initializes the object. 
        """
        assert isinstance(degree, int)
        assert degree > 0
        assert len(supp) == 2
        assert isinstance(supp, list)
        self._degree = degree
        self._measure = pdf
        self._support = supp


    def gram_schmidt(self, grid = 1000):
        """
        Constructs polynomials of order 
        up to "order" that are orthogonal with
        respect to the arbitrary probability 
        measure provided by "ker".
        """
        monos = [lambda x: x**0]
        for i in range(1, self._degree + 1):
            monos = monos + [lambda x: x**i]
        c = np.zeros((self._degree, self._degree))
        dom = np.linspace(self._support[0], self._support[1], grid)
        h = dom[1] - dom[0]
        polys = np.zeros((self._degree+1, grid))
        polys[0,:] = monos[0](dom)
        for i in range(1, self._degree+1):
            for j in range(i):
                c[i-1,j] = np.sum(monos[i](dom) * polys[j,:] * self._measure(dom)) / np.sum(polys[j,:] ** 2 * self._measure(dom))
            polys[i, :] = monos[i](dom) - np.dot(c[i-1,:i], polys[:i,:]) 
        return polys, c


    def recurr_coeffs(self, grid = 1000):
        t = lambda x: x
        poly, c = self.gram_schmidt(grid)
        dom = np.linspace(self._support[0], self._support[1], grid)
        h = dom[1] - dom[0]
        alphas = np.zeros(poly.shape[0])
        betas = np.zeros(poly.shape[0]-1)
        alphas[0] = np.sum(t(dom) * self._measure(dom)) * h
        for i in range(poly.shape[0]-1):
            alphas[i+1] = np.sum( t(dom) * poly[i+1,:]**2 * self._measure(dom)) * h / (np.sum(poly[i+1,:]**2 * self._measure(dom)) * h)
            betas[i] = np.sum(poly[i+1,:]**2 * self._measure(dom)) * h / (np.sum(poly[i,:] ** 2 * self._measure(dom)) * h)
        return alphas, betas


class PolyBasis(object):
    """ construct basis terms matrix """

    _degree = None 
    _dim = None
    _MI_terms = None
    _type = None

    def __init__(self, dim = 1, degree = 1, pol_type = 'H'):
        """
        Ininializes the object 
        """
        assert isinstance(dim, int)
        assert isinstance(degree, int)
        assert dim > 0
        assert degree > 0
        assert pol_type in ['H', 'L', 'Lag'], 'Only Hermite, Legendre and Laguerre polynomials are currently supported ! Choose among ' + str(['H', 'L', 'Lag']) + ' !'
        self._degree = degree 
        self._dim = dim
        self._MI_terms = self.mi_terms(self._dim, self._degree)
        self._type = pol_type


    def __call__(self, XI):
        assert XI.shape[1] == self._dim
        if self._type == 'H':
            H = [Hermite1d(degree = self._degree)(XI[:,i]) for i in range(self._dim)]
            PSI = np.ones((XI.shape[0], self._MI_terms.shape[0]))
            for i in range(self._MI_terms.shape[0]):
                for j in range(self._dim):
                    PSI[:,i] *= H[j][:,self._MI_terms[i,j]]
        elif self._type == 'L':
            H = [Legendre1d(degree = self._degree)(XI[:,i]) for i in range(self._dim)]
            PSI = np.ones((XI.shape[0], self._MI_terms.shape[0]))
            for i in range(self._MI_terms.shape[0]):
                for j in range(self._dim):
                    PSI[:,i] *= H[j][:,self._MI_terms[i,j]]
        else:
            H = [Laguerre1d(degree = self._degree)(XI[:,i]) for i in range(self._dim)]
            PSI = np.ones((XI.shape[0], self._MI_terms.shape[0]))
            for i in range(self._MI_terms.shape[0]):
                for j in range(self._dim):
                    PSI[:,i] *= H[j][:,self._MI_terms[i,j]]
        return PSI


    def mi_terms(self, dim = None, order = None, trunc = 'TD'):
        """ matrix of multi-indices corresponding to basis terms
        Input
        :order: PCE order
        :dim: PCE dimension
        :trunc: truncation type 
            'TD' --> total degree (sum alpha_i <= order)
            'TP' --> tensor product (max alpha_i <= order)
            'LQ' --> l_q truncation (sum alpha_i^q <= order^q)
            'HC' --> hyperbolic cross (prod [alpha_i+1] <= [order + 1])
        """
        assert trunc in ['TD', 'TP'], 'Only total degree (TD) and tensor product (TP) truncation is currently supported ! l_q and hyperbolic cross truncation are underway !'
        
        if dim is None:
            dim = self._dim
        if order is None:
            order = self._degree

        if trunc == 'TD':
            q_num = [int(misc.comb(dim+i-1, i)) for i in range(order+1)]
            mul_ind = np.array(np.zeros(dim, dtype = int), dtype = int)
            mul_ind = np.vstack([mul_ind, np.eye(dim, dtype = int)])
            I = np.eye(dim, dtype = int)
            ind = [1] * dim
            for j in range(1,order):
                ind_new = []
                for i in range(dim):
                    a0 = np.copy(I[int(np.sum(ind[:i])):,:])
                    a0[:,i] += 1
                    mul_ind = np.vstack([mul_ind, a0])
                    ind_new += [a0.shape[0]]
                ind = ind_new
                I = np.copy(mul_ind[np.sum(q_num[:j+1]):,:])
        elif trunc == 'TP':
            x = np.arange(dim+1)
            mul_ind = np.array(list(itls.product(x,x)))[:,::-1]

        return mul_ind
    
    def mi_terms_loc(self, d1, d2, ord):
        assert d1 < d2
        MI = self.mi_terms(d2, ord)
        if d2 == d1 + 1:
            return np.where(MI[:,-1] == 0)[0]
        else:
            inds = np.vstack([np.where(MI[:,d1:] == [0]*(d2-d1))[0], np.where(MI[:,d1:] == [0]*(d2-d1))[1]]).T
            locs = []
            for i in range(MI.shape[0]):
                if len(np.where( inds[:,0] == i)[0]) == d2-d1:
                    locs += [i]
            return locs

