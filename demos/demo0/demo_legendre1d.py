import chaos_basispy as cb
import numpy as np
import scipy.stats as st 
import matplotlib.pyplot as plt

# Creates 1d Legendre polynomials of order up to 50 and evaluates them on Monte Carlo samples

deg = 50

pol1d = cb.Legendre1d(deg)
#xi = np.linspace(-3.5, 3.5, 1000)
xi = st.uniform(loc = -1, scale = 2.).rvs(size = (100000,))
H = pol1d(xi)
fig = plt.figure(figsize = (8,7))
for i in range(16):
	ax = fig.add_subplot(4,4,i+1)
	ax.plot(xi, H[:,i], '.', ms = 1)

fig_norm = plt.figure()
ax_norm = fig_norm.add_subplot(111)
ax_norm.plot(np.std(H, axis = 0), '.', ms = 1.5)
ax_norm.set_xlabel('order n', fontsize = 12)
ax_norm.set_title('Norm', fontsize = 12)
plt.show()
