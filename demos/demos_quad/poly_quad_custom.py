import numpy as np 
import chaos_basispy as cb
import matplotlib.pyplot as plt


ker = lambda x: 1/2.
supp = [-1., 1.]

custom = {'pdf': ker, 'supp': supp}

quad = cb.QuadratureRule('custom_pdf')
quad2 = cb.QuadratureRule('GL')

[x1, w1] = quad.get_rule(2, 10, custom = custom)

[x2, w2] = quad2.get_rule(2, 11, False)

print x1.shape, x2.shape

fig = plt.figure(figsize = (12,5))
ax1 = fig.add_subplot(121)
ax1.plot(x1[:,0], x1[:,1], '.')
ax2 = fig.add_subplot(122)
ax2.plot(x2[:,0], x2[:,1], '.')
plt.show()

