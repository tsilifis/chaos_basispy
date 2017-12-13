import numpy as np 
import scipy.stats as st 
import matplotlib.pyplot as plt 
import chaos_basispy as cb

def f(xi, a, b, c, W):
    assert xi.shape[0] == 10
    assert W.shape[0] == 10
    return a + b * np.dot(W.T, xi) + c * np.dot(xi.reshape(1,xi.shape[0]), np.dot(np.dot(W, W.T) , xi))

dim = 10
np.random.seed(1234)
W = np.random.normal(size = (dim,))
W = W.reshape(W.shape[0], 1) / np.linalg.norm(W)
print W 
a = np.random.normal()
b = np.random.normal()
c = np.random.normal()
print a, b, c

fig = plt.figure(figsize = (9,6))
for l in range(1,4):
	[x, w] = cb.QuadratureRule('CC').get_rule(10, l)
	print '-'*10 +' Level ' + str(l) + ' : ' + str(x.shape[0]) + ' abscissae'

	out = np.zeros(x.shape[0])
	for i in range(x.shape[0]):
		out[i] = f(x[i,:], a, b, c, W)

	pol = cb.PolyBasis(dim, degree = 2, pol_type = 'L')
	P = pol(x)

	coeffs = np.dot(out*w, P)

	xi = st.uniform(loc = -1., scale = 2.).rvs(size = (100,10))
	out_mc = np.zeros(100)
	for i in range(100):
		out_mc[i] = f(xi[i,:], a, b, c, W)

	out_pce = np.dot(pol(xi), coeffs)

	ax1 = fig.add_subplot(2,3,l)
	ax1.plot(x[:,0], x[:,1], 'k.', alpha = 0.8)
	ax1.set_title('Level = '+str(l))
	ax2 = fig.add_subplot(2,3,3+l)
	ax2.plot(out_mc, 'x')
	ax2.plot(out_pce, '.')

plt.show()

