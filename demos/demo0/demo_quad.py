import chaos_basispy as cb
import numpy as np 
import matplotlib.pyplot as plt 

# Initialize a QuadratureRule object for each rule and generate a 2d sparse grid with level 4
rules = ['GH', 'GL', 'NC', 'CC', 'GLag']
dim = 2
lev = 4


fig = plt.figure(figsize = (10,6))
for i in xrange(5):
	quadr = cb.QuadratureRule(rule = rules[i])
	[X, w] = quadr.get_rule(dim, lev)
	ax = fig.add_subplot(2,3,i+1)
	ax.plot(X[:,0], X[:,1], 'o', ms = 3, alpha = 0.8)
	ax.set_title(rules[i], fontsize = 12)

plt.show()
