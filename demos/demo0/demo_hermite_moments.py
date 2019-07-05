import chaos_basispy as cb
import numpy as np
import matplotlib.pyplot as plt


dim = 5
[x, w] = cb.QuadratureRule().get_rule(dim, 6, False)
xi = np.random.normal(size = (100000,5))

pol = cb.PolyBasis(dim, 5)
P_quad = pol(x)
P_MC = pol(xi)
n2_quad = np.dot(w, P_quad**2)
n2_MC = np.std(P_MC, axis = 0)

fig = plt.figure(figsize = (9,7))
ax1 = fig.add_subplot(221)
ax1.plot(x[:,0], x[:,1], 'o', ms = 3, alpha = 0.6)
ax1.set_title('G-H Quadrature ('+str(x.shape[0])+' points)', fontsize = 12)
ax2 = fig.add_subplot(222)
ax2.plot(xi[:,0], xi[:,1], 'o', ms = 2, alpha = 0.3)
ax2.set_title('Monte Carlo ('+r'$10^5$'+' points)', fontsize = 12)
ax3 = fig.add_subplot(223)
ax3.plot(n2_quad, '.', ms = 1.5)
ax3.set_ylim([-0.1,1.2])
ax3.set_xlabel('# of coefficient', fontsize = 12)
ax3.set_ylabel('L2 norm', fontsize = 12)
ax4 = fig.add_subplot(224)
ax4.plot(n2_MC, '.', ms = 1.5)
ax4.set_xlabel('# of coefficient', fontsize = 12)
plt.show()
