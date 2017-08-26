"""
Sparse Grid generator.

Author: Panos Tsilifis
Date: 8/25/2017

"""

__all__ = ['QuadratureRule']


import numpy as np
import math
import itertools as itl
from scipy import misc


class QuadratureRule(object):


	_rule_list = ['GL', 'GH', 'NC', 'CC', 'GLag']
	_rule = None


	@property 
	def rule(self):
		return self._rule

	@rule.setter
	def rule(self, value):
		assert value in self._rule_list, 'Input rule is not (yet) supported. Currently available options are : \n' + str(self._rule_list)
		self._rule = value

	def __init__(self, rule = 'GH'):
		assert rule in self._rule_list, 'Input rule is not (yet) supported. Currently available options are : \n' + str(self._rule_list)
		self._rule = rule

	def mi_terms(self, order, dim):
		""" matrix of basis terms
		Input
		:order: PCE order
		:dim: PCE dimension        
		"""
		q_num = [int(misc.comb(dim+i-1, i)) for i in range(order+1)]
		mul_ind = np.array(np.zeros(dim, dtype = int), dtype = int)
		mul_ind = np.vstack([mul_ind, np.eye(dim, dtype = int)])
		I = np.eye(dim, dtype = int)
		ind = [1] * dim
		terms = []
		for j in range(1,order):
			ind_new = []
			for i in range(dim):
				a0 = np.copy(I[int(np.sum(ind[:i])):,:])
				a0[:,i] += 1
				mul_ind = np.vstack([mul_ind, a0])
				ind_new += [a0.shape[0]]
			ind = ind_new
			I = np.copy(mul_ind[np.sum(q_num[:j+1]):,:])
			terms = [list(mul_ind[i,:]) for i in range(mul_ind.shape[0]) if mul_ind[i,:].min() > 0]
		return np.array(terms)


	def Hermite(self, x,n):
		"""
		Computes the (probabilist's) Hermite polynomials of 
		order up to n at point x. 
		"""
		if n == 1:
			return 1.
		else:
			H = np.zeros(n)
			H[0] = 1.
			H[1] = x
			for i in range(2,n):
				H[i] = x*H[i-1] - (i-1) * H[i-2]
			H = H / [math.sqrt(math.factorial(i)) for i in range(H.shape[0])]
			return H


	def Legendre(self, x,n):
		"""
		Computes the Legendre polynomials of 
		order up to n at point x.
		"""
		if n == 1:
			return 1.
		else:
			H = np.zeros(n)
			H[0] = 1.
			H[1] = x
			for i in range(2, n):
				H[i] = ( (2.*i-1.) * x * H[i-1] - (i - 1.) * H[i-2] ) / i
			H = H / [math.sqrt(1. / (2.*i+1.)) for i in range(H.shape[0])]
			return H

	def Laguerre(self, x,n):
		"""
		Computes the Laguerre polynomials of 
		order up to n at point x.
		"""
		if n == 1:
			return 1.
		else:
			H = np.zeros(n)
			H[0] = 1.
			H[1] = 1. - x
			for i in range(2, n):
			#H[i] = (x - 2*(i-1) - 1.) * H[i-1] - (i-1)**2 * H[i-2]
				H[i] = ((2*(i-1)+1-x)*H[i-1] - (i-1)*H[i-2]) / (i)
			return H


	def GaussHermite(self, n, odd = False):
		assert n > 0
		assert isinstance(n, int)
		if odd:
			n = 2 ** n - 1
		d = np.sqrt(np.arange(n))[1:]
		w = np.zeros(n)
		H = np.diag(d, -1) + np.diag(d, 1)
		[x, v] = np.linalg.eigh(H)
		if odd:
			x[(n-1)/2] = 0.
		for i in range(n):
			w[i] = 1 / np.sum(self.Hermite(x[i], n)**2)
		rule = {'x': x, 'w': w}
		return rule


	def GaussLegendre(self, n, odd = False):
		assert n > 0
		assert isinstance(n, int)
		if odd:
			n = 2 ** n - 1
		d = np.sqrt([i**2/((2.*i+1.)*(2*i-1.)) for i in range(1,n)])
		w = np.zeros(n)
		H = np.diag(d, -1) + np.diag(d, 1)
		[x, v] = np.linalg.eigh(H)
		for i in range(n):
			w[i] = 1 / np.sum(self.Legendre(x[i], n) ** 2)
		if (x.shape[0]-1) % 2 == 0:
			x[(x.shape[0]-1) / 2] = 0.
		rule = {'x': x, 'w': w}
		return rule


	def GaussLaguerre(self, n, odd = False):
		assert n > 0
		assert isinstance(n, int)
		if odd:
			n = 2 ** n - 1
		d = [i for i in range(1,n)]
		alpha = [2*i+1 for i in range(n)]
		w = np.zeros(n)
		H = np.diag(alpha) + np.diag(d, -1) + np.diag(d, 1)
		[x, v] = np.linalg.eigh(H)
		for i in range(n):
			w[i] = 1 / np.sum(self.Laguerre(x[i], n) ** 2)
		rule = {'x': x, 'w': w}
		return rule


	def Trapezoidal(self, l):
		assert l > 0
		assert isinstance(l, int)
		if l == 1:
			x = np.array([0.])
			w = np.array([1.])
		else:
			n = 2 ** (l - 1) + 1
			x = np.array([i*2**(2-l)-1. for i in range(n)])
			w = np.array([2 ** (1.-l) for i in range(n)])
			w[0] = w[0] / 2.
			w[-1] = w[-1] / 2.
		return x, w


	def ClenshawCurtis(self, l):
		assert l > 0
		assert isinstance(l, int)
		if l == 1:
			x = np.array([-np.cos(math.pi / 2.)])
			w = np.array([2.])
		else:
			n = 2 ** (l - 1) + 1
			#n = 2*l - 1
			x = np.array([-np.cos( math.pi * i / (n - 1.)) for i in range(n)])
			w = np.zeros(n)
			w[0] = w[-1] = 1. / (n * (n - 2.))
			for i in range(1, n-1):
				s = [np.cos(2* math.pi * i * (j+1)/(n - 1.)) / (1. - 4*(j+1.)**2) for j in range((n-1)/2)]
				s[-1] = s[-1] / 2.
				w[i] = 2 * (1. + 2. * np.sum(s)) / (n - 1.)
		return x, w / 2.


	def RuleDiff(self, rule1, rule2, nested = 'NN'):
		assert nested in ['NN', 'FN'], 'Rules have to be either Non-nested (NN) or Fully-nested (FN)'
		if nested == 'NN':
			return np.hstack([rule1['x'], rule2['x']]), np.hstack([rule1['w'], -rule2['w']])
		else:
			locs = [i for i in range(len(rule1[0])) for j in range(len(rule2[0])) if rule1[0][i]==rule2[0][j]]
			print locs
			rule1[1][locs] = rule1[1][locs] - rule2[1]
			return rule1[0], rule1[1]


	def get_rule(self, d, l, exp = True):
		assert l > 0
		assert d > 0
		assert isinstance(l, int)
		assert isinstance(d, int)
		#rules = ['CC', 'GL', 'GH', 'NC']
		#assert rule in rules, 'Input rule not supported or non existing.'
		if d == 1:
			if self._rule == 'GH':
				grid = self.GaussHermite(l, exp)
				[x, w] = grid['x'], grid['w']
			elif self._rule == 'GL':
				grid = self.GaussLegendre(l, exp)
				[x, w] = grid['x'], grid['w']
			elif self._rule == 'GLag':
				grid = self.GaussLaguerre(l, exp)
				[x, w] = grid['x'], grid['w']
			elif self._rule == 'NC':
				[x, w] = self.Trapezoidal(l)
			else:
				[x, w] = self.ClenshawCurtis(l)
			return x, w
		else:
			if self._rule in ['CC', 'NC']:
				nested = 'FN'
			elif self._rule in ['GH', 'GL']:
				nested = 'NN'
			d0 = d - 1
			if self._rule == 'GH':
				Q1 = [self.GaussHermite(i, exp) for i in range(1,l+d0)]
				H_rule = self.GaussHermite(1, exp)
				D1 = [[H_rule['x'], H_rule['w']]]
			elif self._rule == 'GL':
				Q1 = [self.GaussLegendre(i, exp) for i in range(1,l+d0)]
				L_rule = self.GaussLegendre(1, exp)
				D1 = [[L_rule['x'], L_rule['w']]]
			elif self._rule == 'GLag':
				Q1 = [self.GaussLaguerre(i, exp) for i in range(1,l+d0)]
				L_rule = self.GaussLaguerre(1, exp)
				D1 = [[L_rule['x'], L_rule['w']]]
			elif self._rule == 'CC':
				Q1 = [self.ClenshawCurtis(i) for i in range(1,l+d0)]
				D1 = [self.ClenshawCurtis(1)]
			elif self._rule == 'NC':
				Q1 = [self.Trapezoidal(i) for i in range(1,l+d0)]
				D1 = [self.Trapezoidal(1)]
			for i in range(1,len(Q1)):
				D1 = D1 + [self.RuleDiff(Q1[i], Q1[i-1])]

			MI = self.mi_terms(l+d0-1, d0)
			THETA = np.zeros(d)
			W = np.zeros(1)
			for i in range(MI.shape[0]):
				mi = MI[i,:]
				k_abs = mi.sum()
				x = [(y,) for y in Q1[l+d0-1-k_abs]['x']]
				w = [(y,) for y in Q1[l+d0-1-k_abs]['w']]
				for j in range(d0):
					x = [u+(v,) for u in x for v in D1[mi[j]-1][0]]
					w = [u+(v,) for u in w for v in D1[mi[j]-1][1]]
				THETA = np.vstack([THETA, np.array(x)])
				W = np.hstack([W, np.prod(np.array(w), axis = 1)])
			if self._rule in ['GH', 'GL'] and exp == True:
				W = np.delete(W, 0,0)
				[theta_uni, ind, inv, c] = np.unique(np.delete(THETA,0,0), True, True, True, axis = 0)
				locs = [j for j in range(c.shape[0]) if c[j] > 1]
				w_uni = W[ind]
				for j in locs:
					loc_j = [k for k in range(inv.shape[0]) if inv[k] == j]
					w_uni[j] = W[loc_j].sum()
				return theta_uni, w_uni
			else:
				return np.delete(THETA, 0, 0), np.delete(W, 0, 0)



