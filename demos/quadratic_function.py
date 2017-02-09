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