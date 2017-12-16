"""

"""


__all__ = ['PolyChaos', 'BasisAdapt']


try:
    import numpy as np
except ImportError:
    print ('numpy needs to be installed')

try:
    import scipy.stats as st 
    from scipy.linalg import polar
except ImportError:
    print ('scipy needs to be installed')

from _poly_basis import Hermite1d, PolyBasis


class PolyChaos(object):
    """
    Class for computing Polynomial Chaos coefficients
    and for generating corresponding samples
    """

    _name = None

    def __init__(self, name = 'PC expansion'):
        """
        Initializes the object
        """
        self._name = name


    def comp_coeffs(self, xi, y, w, degree = 1, poly_type = 'H'):
        """
        Computes Chaos coefficients by performing numerical integration 
        on the expression c_{alpha} = E[y psi_{alpha}(xi)] based on quadrature rules.
        Input
        :xi: array (N,d) where N the number of quadrature points and 
                               and d the dimension of the input variables.
        :y: array (N,) the model output values (should be of equal length with xi.shape[0]).
        :w: array (N,) the quadrature weights.
        :degree: int
                The highest order of the polynomials in the expansion.
        :poly_type: str 
               The type of polynomials in the expansion (currently supports only Hermite and Legendre)

        [1] Tipireddy, R. and Ghanem, R., 2014. Basis adaptation in homogeneous chaos spaces. 
        Journal of Computational Physics, 259, pp.304-317.
        """
        assert xi.shape[0] == y.shape[0]
        assert xi.shape[0] == w.shape[0]
        assert isinstance(degree, int)
        d = xi.shape[1]
        if poly_type not in ['H', 'L']:
            raise ValueError('Unknown type of polynomials !')
        pol = PolyBasis(d, degree, poly_type)
        return np.dot(y * w, pol(xi))


    def sample(self, coeffs, dim, degree, N_mc = 1e+3, poly_type = 'H'):
        """  
        Generates samples from a Polynomial Chaos expansion whose coefficients and type 
        of polynomials are provided.
        """
        assert N_mc > 0
        assert isinstance(N_mc, int)

        pol = PolyBasis(dim , degree)
        assert pol._MI_terms.shape[0] == coeffs.shape[0]

        if poly_type == 'H':
            X = st.norm.rvs(size = (N_mc, dim))
        elif poly_type == 'L':
            X = st.uniform(loc = -1., scale = 2.).rvs(size = (N_mc, dim))
        else:
            raise ValueError('Unknown type of polynomials !')
        return np.dot(pol(X), coeffs)


    def sample_many_QoI(self, coeffs, orders, dim, MC1 = 1e5, X1 = [213], ptype = 'Hermite', seed = 1):
        """
        Generates Monte Carlo samples for multiple QoIs for which the PCE coefficients are provided. Allows 
        PCEs to have different highest order.

        Input: 
                coeffs: list that contains numpy arrays with shape (NC_i,) where NC_i the number of coeffs for the ith QoI.
                orders: list that contains integers specifying the order of the ith QoI.
                dim : int, the input dimensionality
                The rest as in generate_samples_from_pce
        """
        assert len(coeffs) == len(orders)
        m = len(coeffs)
        samples = np.zeros((MC1,m))
        for i in range(m):
            [samples[:,i], xi] = self.generate_samples_from_pce(self, coeffs[i],dim=dim,order=orders[i],MC1=MC1,X1=X1,ptype=ptype, seed = seed)
        return np.vstack([xi, samples])



    def l2_error(self, q_1, q_2, d1, d2, order):
        """ l2-norm relative error function
        Input
        :q_1: vector of coefficients corresponding to eta_{i-1}
        :q_2: vector of coefficients corresonding to eta_{i-2}
        Return
        :e: relative l2-error ||q_1 - q_2|| / ||q_2||   

        """
        import math
        assert np.shape(q_2)[0] == math.factorial(d2+order) / (math.factorial(order) * math.factorial(d2))
        assert np.shape(q_1)[0] < np.shape(q_2)[0]
        Q1 = np.zeros(q_2.shape[0])
        #mi = PolyBasis().mi_terms(order, d2)
        Q1[PolyBasis().mi_terms_loc(d1, d2, order)] = q_1
        return (np.linalg.norm((Q1 - q_2),2) / np.linalg.norm(q_2,2))


    def sobol_ind(self, idx, od, dim, C):      
        """ Computate PC-based 
        Sobol sensitivity
        indices 
        
        Input
        :idx: index of xi for which to compute
              Sobol's sensitivity indices 
              (start from 0)
        :od: order of highest polynomial of PCE
        :dim: dimensionality
        :C: matrix of PCE coefficients
            corresponding to the multi-
            index matrix
        
        Return 
        :su: Sobol PC-based sensitivity indices 
             corresponding to xi's at indices idx

        Ref: 
            Sudret (2006) 'Global sensitivity analysis
            using polynomial chaos expansion'  
        """
        assert isinstance(idx, list)
        if len(C.shape) == 1:
            C = C.reshape(C.shape[0],1)

        MI = PolyBasis().mi_terms(dim, od)
        assert MI.shape[0] == C.shape[0]
        loc = range(MI.shape[0])

        for i in np.arange(dim):
            ai = None
            if i in idx:
                ai = np.where(MI[:,i] > 0)
            else:
                ai = np.where(MI[:,i] == 0)
            loc = np.intersect1d(loc,ai)
            
        Cx = C[loc,:]  
        su = np.zeros(Cx.shape[1])
        
        for m in range(np.shape(Cx)[1]):
            su[m] = np.sum(Cx[:,m]**2) / np.sum(C[1:,m]**2)
        
        return su


class BasisAdapt(PolyChaos):
    """ Basis adaptation """


    def __init__(self, name = 'Adapted PC expansion'):
        """
        Ininialize the object
        """

        super(BasisAdapt, self).__init__(name = name)


    def gauss_adaptation(self, coeffs, dim, method = 0):
        """
        Computes the rotation matrix A that corresponds to the Gaussian adaptation [1]. 
        Input
        :coeffs : array 
            The coefficients of the first order Hermite polynomials (linear part)
            that consist the first row of A.
        :dim: int 
            Dimension of the original input germs xi.
        :method: int {0,1,2,3}
            The method to compute the isometry. 
                0 (default) : Via orthogonal decomposition of a * a.T
                1 : Via a Gram-Schmidt procedure on the matrix C with coeffs (normalized) at its first row,
                        ones along its diagonal, zeros anywhere else.
                2 : Via orthogonal decomposition of the Householder matrix H = I - 2 u u' / ||u|| ** 2
                3 : Via a Gram-Schmidt procedure on the matrix C with coeffs (normalized) at its first row,
                        and the ith largest coeff on the column corresponding to its variable xi
                        on the (i+1) row, zero elsewhere
                4 : Via a Gram-Schmidt procedure on the matrix C with N-QoI coeffs (normalized) at its 
                        first N rows, ones along its diagonal, zeros anywhere else.
                5 : Via a polar decomposition of matrix C with coeffs (normalized) at its first row,
                        and the ith largest coeff on the column corresponding to its variable xi
                        on the (i+1) row, zero elsewhere 
        [1] Tipireddy, R. and Ghanem, R., 2014. Basis adaptation in homogeneous chaos spaces. 
        Journal of Computational Physics, 259, pp.304-317.
        """

        assert coeffs.shape[0] == dim
        if method != 4:
            a = coeffs.reshape((coeffs.shape[0], 1))
        elif method == 4:
            a = coeffs.reshape((coeffs.shape[0], coeffs.shape[1]))
        else:
            print ('Shape of first order coefficients unrecognized')

        #H = np.eye(dim) - 2 * np.dot(u, u.T) / np.linalg.norm(u) ** 2 # Householder matrix
        if method == 0:
            C = np.dot(a, a.T)
            [vals, vecs] = np.linalg.eigh(C)
            return vecs[:,::-1].T
        elif method == 1: # Gram-Schmidt
            C = np.eye(a.shape[0])
            C[0,:] = a[:,0]
            [q,r] = np.linalg.qr(C.T)
            return q.T
        elif method == 2: # Householder matrix
            H = np.eye(dim) - 2 * np.dot(a, a.T) / np.linalg.norm(a) ** 2 # Householder matrix
            [vals, vecs] = np.linalg.eigh(H)
            return np.real(vecs).T
        elif method == 3:
            c3 = np.argsort(np.abs(coeffs))[::-1]
            C = np.zeros((dim, dim))
            C[:,0] = a[:,0]
            j = 0        
            for i in range(0, np.size(coeffs)-1):
                C[c3[j], i+1] = coeffs[c3[j]]
                j += 1
            [q,r] = np.linalg.qr(C)
            q=np.mat(q)
            return q.T
        elif method == 4: # Multi-QoI with Gram-Schmidt
            C = np.eye(a.shape[0])
            for k in np.arange(coeffs.shape[1]):
                C[k,:] = a[:,k]    
                [q,r] = np.linalg.qr(C.T)
            return q.T
        elif method == 5: # Polar decomposition
            c3 = np.argsort(np.abs(coeffs))[::-1]
            C = np.zeros((dim, dim))
            C[:,0] = a[:,0]
            j = 0        
            for i in range(0, np.size(coeffs)-1):
                C[c3[j], i+1] = coeffs[c3[j]]
                j += 1
            [u,p] = polar(C)
            return u.T, np.linalg.det(C)            
        else:
            raise ValueError('Method parameter must be in {0,1,2,3,4,5}')

    def quadratic_adaptation(self, coeffs, dim):
        """
        Computes what is referred to as Quadratic Adaptation in [1], using chaos coefficients of the 
        2nd ordre Hermite polynomials.
        Input
        :coeffs: array
                The coefficients of the second order Hermite polynomials (quadratic part)
        :dim: int
                Dimension of the input germs xi.

        [1] Tipireddy, R. and Ghanem, R., 2014. Basis adaptation in homogeneous chaos spaces. 
        Journal of Computational Physics, 259, pp.304-317.  
        """

        assert coeffs.shape[0] == dim * (dim + 1) / 2
        Q_upper = np.zeros((dim, dim))
        for l in range(dim):
            for k in range(dim-l):
                if k == 0:
                    Q_upper[l,l+k] = coeffs[np.sum(np.arange(dim)[dim-l:]) + k] / np.sqrt(2.)
                else:
                    Q_upper[l,l+k] = coeffs[np.sum(np.arange(dim)[dim-l:]) + k] 
        S = Q_upper + np.transpose(Q_upper) - np.diag(np.diag(Q_upper))
        [vals, vecs] = np.linalg.eigh(S, UPLO = 'U')
        return vecs.T

    def eta_to_xi_mapping(self, eta, A, zeta = None):
        """
        Maps points from the (lower dimensional) eta-space to the xi space. The dimensionality of eta
        must be smaller than or equal to that of xi (shape of A) and the rotation matrix must be a 
        unitaty matrix.
        Parameters ---- > eta : array (N, d0) where N the number of points, d0 the dimensionality of eta
                                The set of N (quadrature) points on a d0-dimensional space.
                            A : array (d, d)
                                Isometry which serves as the rotation matrix (xi = A' [eta, zeta]) 
                         zeta : array (N, d-d0) 
                                Optional set of N points on a (d-d0)-dimensional space that will 
                                augment eta to the d-dimensional space. Zero values will be used at the 
                                default choice when not user-specified.

        [1] Tipireddy, R. and Ghanem, R., 2014. Basis adaptation in homogeneous chaos spaces. 
        Journal of Computational Physics, 259, pp.304-317.
        """

        assert A.shape[0] == A.shape[1]

        d0 = eta.shape[1]
        d = A.shape[0]
        N = eta.shape[0]
        if zeta == None:
            zeta = np.zeros((N, d-d0))
        else:
            assert eta.shape[0] == zeta.shape[0]
            assert eta.shape[1] + zeta.shape[1] == A.shape[0]
        eta_full = np.hstack([eta, zeta]) # Augment eta points with zeta
        return np.dot(A.T, eta_full.T).T


    def transform_coeffs(self, coeffs, deg, iso, eta_dim = 1, method = 0, num_MC = 10000, xi = None, w = None):
        """
        Given the coefficients of a (low dimensional) chaos expansion with respect to the 
        eta basis, where eta = A * xi, transforms the coefficients to those that correspond 
        to an expansion with respect to the original xi basis.
        Input
        :coeffs: array(n,) 
                The coefficients of the PCE wrt the eta basis
        :iso: array(d,d) where d is the dimensionality of xi.
               The (unitary) matrix that relates eta = A * xi.
        :eta_dim: int
                The number of eta components used in the current PCE (eta_dim <= d).
        :num_MC: int
                The number of Monte Carlo samples to be used for estimating the 
                inner products <psi_{beta}(A * xi), psi_{alpha}(xi)>

        [1] Tipireddy, R. and Ghanem, R., 2014. Basis adaptation in homogeneous chaos spaces. 
        Journal of Computational Physics, 259, pp.304-317.
        """

        assert isinstance(eta_dim, int)
        assert isinstance(deg, int)
        assert eta_dim <= iso.shape[0]
        pol_eta = PolyBasis(eta_dim, deg)
        assert len(pol_eta._terms) == coeffs.shape[0]
        if method == 0:
            xi = st.norm.rvs(size = (num_MC, iso.shape[1]))
            eta = np.dot(iso, xi.T)[:eta_dim, :].T
            pol = PolyBasis(iso.shape[0], deg)
            return np.dot(coeffs, np.dot(pol_eta(eta).T, pol(xi))) / num_MC
        elif method == 1:
            eta = np.dot(iso, xi.T)[:eta_dim, :].T
            pol = PolyBasis(iso.shape[0], deg)
            return np.dot(coeffs, np.dot(pol_eta(eta).T * w, pol(xi)))
        else:
            raise ValueError('For integration with Monte Carlo or quadrature rule choose 0 or 1 respectively.')



