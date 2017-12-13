import numpy as np 
import chaos_toolbox as ct
import matplotlib.pyplot as plt

# Uniform(-1, 1)
ker1 = lambda x: 1/2.
supp1 = [-1., 1.]
custom1 = {'pdf': ker1, 'supp': supp1}

# Parabolic pdf on (-1, 1)
ker2 = lambda x: x**2 * 3/2.
custom2 = {'pdf': ker2, 'supp': supp1}

# Cubic pdf on (-1, 1)
ker3 = lambda x: x**3 / 2. + 1/2.
custom3 = {'pdf': ker3, 'supp': supp1}

# Uniform(-1.5, 1.5)
ker4 = lambda x: 5 * x ** 4 / 2.
custom4 = {'pdf': ker4, 'supp': supp1}

quad = ct.util.QuadratureRule('custom_pdf')
[x1, w1] = quad.get_rule(2, 10, custom = custom1)
[x2, w2] = quad.get_rule(2, 10, custom = custom2)
[x3, w3] = quad.get_rule(2, 10, custom = custom3)
[x4, w4] = quad.get_rule(2, 10, custom = custom4)


#print x1
#print x2
#print w
#print w.sum()

#plt.plot(x, '.')
#plt.plot(x2, '*')
#plt.show()



#locs = [i for i in range(x.shape[0]) if np.abs(x[i,:]).max() > 1.]
#locs_c = list(set(range(x.shape[0])) - set(locs))
#y = x[locs_c,:]

print w1.sum(), w2.sum(), w3.sum(), w4.sum()

fig = plt.figure(figsize = (10,12))
ax1 = fig.add_subplot(221)
ax1.plot(x1[:,0], x1[:,1], '.')
ax1.set_title('Uniform', fontsize = 12)
ax2 = fig.add_subplot(222)
ax2.plot(x2[:,0], x2[:,1], '.')
ax2.set_title('Parabolic', fontsize = 12)
ax3 = fig.add_subplot(223)
ax3.plot(x3[:,0], x3[:,1], '.')
ax3.set_title('Qubic', fontsize = 12)
ax4 = fig.add_subplot(224)
ax4.plot(x4[:,0], x4[:,1], '.')
ax4.set_title('Quartic', fontsize = 12)
plt.show()

