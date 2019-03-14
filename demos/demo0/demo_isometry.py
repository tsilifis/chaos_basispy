"""
Visualize isometry constructed 
using different orthogonalization methods

> d = 15
> random coefficients qe ~ N(0,1)

"""
import numpy as np 
import matplotlib.pyplot as plt
import chaos_basispy as cb
import matplotlib
matplotlib.rcParams.update({'font.size': 7})


def visualize_isometry(d):
    d
    qe = np.random.normal(size=(d,))
    
    A0 = cb.BasisAdapt().gauss_adaptation(qe, d, method = 0)
    A1 = cb.BasisAdapt().gauss_adaptation(qe, d, method = 1)
    A2 = cb.BasisAdapt().gauss_adaptation(qe, d, method = 2)
    
    fig = plt.figure() 
    ax = fig.add_subplot(1,3,1)
    ax.imshow(A0, interpolation='none')
    plt.title(r'$method=0$',fontsize = 9)
    
    ax = fig.add_subplot(1,3,2)
    ax.imshow(A1, interpolation='none')
    plt.title(r'$method=1$',fontsize = 9)
    
    ax = fig.add_subplot(1,3,3)
    ax.imshow(A2, interpolation='none')
    plt.title(r'$method=2$', fontsize = 9)
    
    plt.savefig('isometry.eps', dpi=1200)

if __name__ == '__main__':
    visualize_isometry(d = 15)
    plt.show()
