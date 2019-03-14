"""
Construct a PCE for the solution of 
the diffusion PDE.

Author:
    Panagiotis Tsilifis

Date:
    1/10/2018

"""

import numpy as np
import fipy as fp
import chaos_basispy as cb
import matplotlib.pyplot as plt


# Make the source 

def forward(xs):

    nx = 21
    ny = nx
    dx = 1./51 
    dy = dx
    rho = 0.05
    q0 = 1. / (np.pi * rho ** 2)
    T = 0.3
    mesh = fp.Grid2D(dx=dx, dy=dy, nx=nx, ny=ny)


    time = fp.Variable()
    sourceTerm_1 = fp.CellVariable(name = "Source term", mesh=mesh, value = 0.)

    for i in range(sourceTerm_1().shape[0]):
        sourceTerm_1()[i] = q0 * np.exp( - ((mesh.cellCenters[0]()[i] - xs[0]) ** 2 
                                  + (mesh.cellCenters[1]()[i] - xs[1]) ** 2 ) / (2 * rho **2)) * (time() < T)

    # The equation
    eq = fp.TransientTerm() == fp.DiffusionTerm(coeff=1.) + sourceTerm_1# + sourceTerm_2

    # The solution variable
    phi = fp.CellVariable(name = "Concentration", mesh=mesh, value=0.)


    #if __name__ == '__main__':
    #    viewer = fp.Viewer(vars=phi, datamin=0., datamax=3.)
    #    viewer.plot()

    x = np.arange(0,nx)/nx
    y = x

    data = []
    dt = 0.005
    steps = 60
    for step in range(steps):
        time.setValue(time() + dt)
        eq.solve(var=phi, dt=dt)
        #if __name__ == '__main__':
        #    viewer.plot()
        #if step == 14 or step == 29 or step == 44 or  step == 59:
        #    dl = phi()[0]
        #    dr = phi()[nx-1]
        #    ul = phi()[nx**2 - nx]
        #    ur = phi()[nx**2 - 1]
        #    print phi().shape
        #    data = np.hstack([data, np.array([dl, dr, ul, ur])])

    #if __name__ == '__main__':
    #    raw_input("Transient diffusion with source term. Press <return> to proceed")

    return phi().reshape(nx, nx)



#plt.plot(x_quad[:,0], x_quad[:,1], '*')
#plt.show()
max_lev = 8
for l in range(2,max_lev+1):
    print '-'*5 + 'Running model at quadrature rule level '+ str(l) + ' points.' + '-'*5
    [x_quad, w_quad] = cb.QuadratureRule('CC').get_rule(2, l)

    u_quad = np.zeros((51, 51, x_quad.shape[0]))

    for i in range(x_quad.shape[0]):
        print 'Run ' + str(i+1)
        u_quad[:,:,i] = forward((x_quad[i,:] + np.ones(2)) / 2)

    np.save('u_quad_lev_'+str(l)+'.npy', u_quad)