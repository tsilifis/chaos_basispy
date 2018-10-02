import numpy as np 
import chaos_basispy as cb
import matplotlib.pyplot as plt 

quad = cb.QuadratureRule('CC')
chaos = cb.PolyChaos()

[x2, w2] = quad.get_rule(2, 2)
[x3, w3] = quad.get_rule(2, 3)
[x4, w4] = quad.get_rule(2, 4)
[x5, w5] = quad.get_rule(2, 5)
[x6, w6] = quad.get_rule(2, 6)
[x7, w7] = quad.get_rule(2, 7)
[x8, w8] = quad.get_rule(2, 8)

x_quads = [x2, x3, x4, x5, x6, x7, x8]
w_quads = [w2, w3, w4, w5, w6, w7, w8]

u2 = np.load('u_quad_lev_2.npy')[25,25,:]
u3 = np.load('u_quad_lev_3.npy')[25,25,:]
u4 = np.load('u_quad_lev_4.npy')[25,25,:]
u5 = np.load('u_quad_lev_5.npy')[25,25,:]
u6 = np.load('u_quad_lev_6.npy')[25,25,:]
u7 = np.load('u_quad_lev_7.npy')[25,25,:]
u8 = np.load('u_quad_lev_8.npy')[25,25,:]

u_quads = [u2, u3, u4, u5, u6, u7, u8]

nx = u2.shape[0]
basis = cb.PolyBasis()

coeffs2 = np.zeros((7, basis.mi_terms(2, 2).shape[0]))
coeffs3 = np.zeros((7, basis.mi_terms(2, 3).shape[0]))
coeffs4 = np.zeros((7, basis.mi_terms(2, 4).shape[0]))
coeffs5 = np.zeros((7, basis.mi_terms(2, 5).shape[0]))
coeffs6 = np.zeros((7, basis.mi_terms(2, 6).shape[0]))
coeffs7 = np.zeros((7, basis.mi_terms(2, 7).shape[0]))
coeffs8 = np.zeros((7, basis.mi_terms(2, 8).shape[0]))
coeffs9 = np.zeros((7, basis.mi_terms(2, 9).shape[0]))
coeffs10 = np.zeros((7, basis.mi_terms(2, 10).shape[0]))

x = np.linspace(0., 1, 51)
y = x
xx, yy = np.meshgrid(x, y)
XI = np.vstack([xx.flatten(), yy.flatten()]).T
print XI.shape
u_chaos_ord2 = np.zeros((x.shape[0], y.shape[0], 7))
for i in range(7):
	coeffs2[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 2, 'L')
	u_chaos_ord2[:,:,i] = np.dot(cb.PolyBasis(2, 2)(XI), coeffs2[i,:]).reshape(51, 51)
	coeffs3[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 3, 'L')
	coeffs4[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 4, 'L')
	coeffs5[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 5, 'L')
	coeffs6[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 6, 'L')
	coeffs7[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 7, 'L')
	coeffs8[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 8, 'L')
	coeffs9[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 9, 'L')
	coeffs10[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 10, 'L')

fig = plt.figure(figsize = (10,7))
ax1 = fig.add_subplot(231)
ax1.contourf(x, y, u_chaos_ord2[:,:,0], 50)
ax2 = fig.add_subplot(232)
ax2.contourf(x, y, u_chaos_ord2[:,:,1], 50)
ax3 = fig.add_subplot(233)
ax3.contourf(x, y, u_chaos_ord2[:,:,2], 50)
ax4 = fig.add_subplot(234)
ax4.contourf(x, y, u_chaos_ord2[:,:,3], 50)
ax5 = fig.add_subplot(235)
ax5.contourf(x, y, u_chaos_ord2[:,:,4], 50)
ax6 = fig.add_subplot(236)
ax6.contourf(x, y, u_chaos_ord2[:,:,5], 50)
#ax7 = fig.add_subplot(237)
#ax7.contourf(x, y, u_chaos_ord2[:,:,6], 50)
plt.show()

plt.plot(coeffs2[0,:], '-x')
plt.plot(coeffs2[1,:], '-s')
plt.plot(coeffs2[2,:], '-+')
plt.plot(coeffs2[3,:], '-d')
plt.plot(coeffs2[4,:], '-*')
plt.plot(coeffs2[5,:], '-^')
plt.plot(coeffs2[6,:], '--.')
plt.show()
