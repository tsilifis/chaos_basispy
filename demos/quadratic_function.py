import numpy as np
import scipy.stats as st 
import matplotlib.pyplot as plt 
import chaospy as cp
import chaos_basispy as cb

dim = 10
np.random.seed(1234)
W = np.random.normal(size = (dim,))
W = W.reshape(W.shape[0],1) / np.linalg.norm(W)
print W
a = np.random.normal()
b = np.random.normal()
c = np.random.normal()
print a,b,c

xi = st.uniform.rvs(loc= -1., scale = 2., size = (10,1000))
#xi = st.norm.rvs(size = (10, 10000))
def f(xi, a, b, c, W):
    assert xi.shape[0] == 10
    assert W.shape[0] == 10
    return a + b * np.dot(W.T, xi) + c * np.dot(xi.reshape(1,xi.shape[0]), np.dot(np.dot(W, W.T) , xi))


out = np.zeros(1000)
for i in range(1000):
    out[i] = f(xi[:,i], a, b, c, W)