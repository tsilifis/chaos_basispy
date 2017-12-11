"""
Demo 1: Illustrates the basis adaptation technique with a toy problem.
The objective function is a 10-dimensional quadratic polynomial with 
known active subspace. 

This example does nor really make use of the module's capabilities but rather
shows the idea behind Basis Adaptation. 

Author : Panagiotis Tsilifis
Date : 2/8/2017
"""

import numpy as np
import scipy.stats as st 
import matplotlib.pyplot as plt 
import chaos_basispy as cb

#============ The quadratic function ==============
def f(xi, a, b, c, W):
    assert xi.shape[0] == 10
    assert W.shape[0] == 10
    return a + b * np.dot(W.T, xi) + c * np.dot(xi.reshape(1,xi.shape[0]), np.dot(np.dot(W, W.T) , xi))


# Fix random seed and generate the projection vector.
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


# Plot the pdf of the projected input eta. 
out = np.zeros(1000)
for i in range(1000):
    out[i] = f(xi[:,i], a, b, c, W)

z = np.dot(W.T, xi)

eta = np.dot(W.flatten(), st.uniform.rvs(loc = -1., scale = 2., size = (dim, 100000)))

plt.style.use('ggplot')

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.hist(eta, bins = 100, normed = True, histtype = 'stepfilled')
ax1.set_xlabel('Value', fontsize = 15)
ax1.set_ylabel('pdf', fontsize = 15)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)  
ax1.get_xaxis().tick_bottom()  
ax1.get_yaxis().tick_left() 
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(15)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)


# Plot the quadrature points on [-1,1] and the images in the projected space.

zeta = np.load('grid_points_CC_dim1_level_5.npy')
weights = np.load('weights_CC_dim1_level_5.npy')

new_eta1 = cb.inverse_transform_sampling(eta, eval_points = (zeta[:,0] + 1.) /2., n_bins=50)


fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(zeta[:,0], zeta[:,0], '.', ms = 10, mew = 0, label = 'CC rule')
ax1.plot(zeta[:,0], new_eta1, '*', ms = 10, mew = 0, label = 'Transformed quadr points')
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)  
ax1.get_xaxis().tick_bottom()  
ax1.get_yaxis().tick_left() 
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(15)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)
plt.legend(loc = 2)


## Construct the 1D PCE and a 10D for comparison. 

new_quad_points = np.zeros((dim,new_eta1.shape[0]))
for i in range(new_eta1.shape[0]):
    new_quad_points[:,i] = W[:,0] * new_eta1[i]


quad_out = np.zeros(new_eta1.shape[0])
for i in range(new_eta1.shape[0]):
    quad_out[i] = f(new_quad_points[:,i], a, b, c, W)

#import orthpol

rvs = [st.uniform(loc = -1., scale = 2.)]
#pol = orthpol.ProductBasis(rvs, degree = 20)
pol = cb.PolyBasis(1, 20, 'L')

print weights.sum(), weights.shape

coeffs = np.zeros(pol.mi_terms(1, 20).shape[0])
P = pol(zeta)

for i in range(weights.shape[0]):
    coeffs = coeffs + quad_out[i] * weights[i] * P[i,:] #/ 2.


zeta_sam = rvs[0].rvs(size = (1000,1))
eta_sam = cb.inverse_transform_sampling(eta, eval_points = (zeta_sam[:,0] + 1.) /2., n_bins=50)
samples = np.dot(pol(zeta_sam), coeffs)

new_xi = np.zeros((10, eta_sam.shape[0]))
for i in range(eta_sam.shape[0]):
    new_xi[:,i] = W[:,0] * eta_sam[i]


quad_points = np.load('grid_points_CC_dim10_level_2.npy')
wghts10 = np.load('weights_CC_dim10_level_2.npy')


quad_out10 = np.zeros(quad_points.shape[0])
for i in range(quad_out10.shape[0]):
    quad_out10[i] = f(quad_points[i,:], a, b, c, W)

#pol10d = orthpol.ProductBasis([st.uniform(loc = -1, scale = 2)]*dim, degree = 2)
pol10d = cb.PolyBasis(dim, 2, 'L')
P10 = pol10d(quad_points)
coeffs10 = np.zeros(pol10d.mi_terms(dim, 2).shape[0])

for i in range(wghts10.shape[0]):
    coeffs10 = coeffs10 + quad_out10[i] * wghts10[i] * P10[i,:]

out10 = np.dot(pol10d(new_xi.T), coeffs10)


et = np.linspace(-2., 2., 1000)

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(et, a + b*et + c * et**2, 'g-', label = r'$g(\eta)$' )
ax1.plot(eta_sam, out10 , 'ko', ms = 3, label = '10d-2ord PCE')
ax1.plot(eta_sam, samples, 'x' , alpha = 0.8, label = '1d-20rd PCE')
ax1.set_xlabel(r'$\eta$', fontsize = 16)
ax1.set_ylabel(r'$g(\eta)$', fontsize = 16)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)  
ax1.get_xaxis().tick_bottom()  
ax1.get_yaxis().tick_left() 
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(15)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)
plt.legend(loc = 2)
#plt.show()

g_range = np.linspace(0., 6., 1001)
ker1 = st.gaussian_kde(out10)
ker2 = st.gaussian_kde(samples)
fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.plot(g_range, ker1(g_range), '-', label = r'$g(\eta)$')
ax1.plot(g_range, ker2(g_range), '-', label = '1d-20ord PCE')
ax1.set_xlabel('Value', fontsize = 15)
ax1.set_ylabel('pdf', fontsize = 15)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)  
ax1.get_xaxis().tick_bottom()  
ax1.get_yaxis().tick_left() 
for tick in ax1.xaxis.get_major_ticks():
    tick.label.set_fontsize(15)
for tick in ax1.yaxis.get_major_ticks():
    tick.label.set_fontsize(15)
plt.legend()
plt.show()

