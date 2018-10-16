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

loc_x, loc_y = 0, -1 # The spatial location of the QoI

u2 = np.load('u_quad_lev_2.npy')[loc_x,loc_y,:]
u3 = np.load('u_quad_lev_3.npy')[loc_x,loc_y,:]
u4 = np.load('u_quad_lev_4.npy')[loc_x,loc_y,:]
u5 = np.load('u_quad_lev_5.npy')[loc_x,loc_y,:]
u6 = np.load('u_quad_lev_6.npy')[loc_x,loc_y,:]
u7 = np.load('u_quad_lev_7.npy')[loc_x,loc_y,:]
u8 = np.load('u_quad_lev_8.npy')[loc_x,loc_y,:]

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
u_chaos_ord3 = np.zeros((x.shape[0], y.shape[0], 7))
u_chaos_ord4 = np.zeros((x.shape[0], y.shape[0], 7))
u_chaos_ord5 = np.zeros((x.shape[0], y.shape[0], 7))
u_chaos_ord6 = np.zeros((x.shape[0], y.shape[0], 7))

for i in range(7):
	coeffs2[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 2, 'L')
	u_chaos_ord2[:,:,i] = np.dot(cb.PolyBasis(2, 2)(XI), coeffs2[i,:]).reshape(51, 51)
	coeffs3[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 3, 'L')
	u_chaos_ord3[:,:,i] = np.dot(cb.PolyBasis(2, 3)(XI), coeffs3[i,:]).reshape(51, 51)
	coeffs4[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 4, 'L')
	u_chaos_ord4[:,:,i] = np.dot(cb.PolyBasis(2, 4)(XI), coeffs4[i,:]).reshape(51, 51)
	coeffs5[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 5, 'L')
	u_chaos_ord5[:,:,i] = np.dot(cb.PolyBasis(2, 5)(XI), coeffs5[i,:]).reshape(51, 51)
	coeffs6[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 6, 'L')
	u_chaos_ord6[:,:,i] = np.dot(cb.PolyBasis(2, 6)(XI), coeffs6[i,:]).reshape(51, 51)
	coeffs7[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 7, 'L')
	coeffs8[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 8, 'L')
	coeffs9[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 9, 'L')
	coeffs10[i, :] = chaos.comp_coeffs(x_quads[i], u_quads[i], w_quads[i], 10, 'L')


coeffs_all = [coeffs2, coeffs3, coeffs4, coeffs5, coeffs6, coeffs7, coeffs8, coeffs9, coeffs10]

fig, axes = plt.subplots(nrows = 2, ncols = 3, figsize = (12, 6))
i = 0
for ax in axes.flat:
	im = ax.contourf(x, y, u_chaos_ord6[:,:,i+1], 50)
	i = i + 1

fig.subplots_adjust(right = 0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(im, cax = cbar_ax)
plt.show()

plt.plot(coeffs4[0,:], '-x')
plt.plot(coeffs4[1,:], '-s')
plt.plot(coeffs4[2,:], '-+')
plt.plot(coeffs4[3,:], '-d')
plt.plot(coeffs4[4,:], '-*')
plt.plot(coeffs4[5,:], '-^')
plt.plot(coeffs4[6,:], '--.')
plt.show()

### --- Estimate truncation error ---

err = np.zeros((7, 9))
for i in range(7):
	for j in range(9):
		err[i,j] = (( u_quads[-1] - np.dot(cb.PolyBasis(2, j+2)(x_quads[-1]), coeffs_all[j][i,:]) )**2 * w_quads[-1]).sum() / (u_quads[-1]**2 * w_quads[-1]).sum()

lev = [2, 3, 4, 5, 6, 7, 8]
Q = [2, 3, 4, 5, 6, 7, 8, 9, 10]

#fig = plt.figure()
#ax = fig.add_subplot(111)
plt.contour(Q, lev, err, 50)
plt.colorbar()
plt.show()

