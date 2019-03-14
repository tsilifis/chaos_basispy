# chaos_basispy
A Polynomial Chaos Basis Reduction module
=======================================

Author:       Panagiotis Tsilifis

Affiliation:  CSQI, Institure of Mathematics, 
              École Polytechnique Fédérale de Lausanne, Lausanne, CH-1015, Switzerland

email address: panagiotis.tsilifis@epfl.ch

Description
-----------

The package 'chaos_basispy' defines the 'chaos_basispy' module which attempts to 
collect the recently developed Basis Adaptation (BA) techniques that have been developed 
particularly for Polynomial Chaos expansions (PCE). 

The core idea behind Basis Adaptation is to apply a transformation on the input variables
in a way such that the Quantity of Interest (QoI) which we are approximating via a PCE can
be expressed with respect to a reduced basis. Provided that the input variables are Gaussian
and the transformation is through an isometry (unitary matrix), the new variables preserve 
"Gaussianity" and a new PCE can be constructed with respect to them. This makes the BA 
framework particularly attractive on Hermite PCEs and in fact inapplicable (yet) to 
generalized PCEs. An expeption is when the new variables is 1-dimensional which can be mapped
to a uniform r.v. through its own cdf and therefore the Legendre adapted PCE can be constructed.
Due to the above, this module considers only Hermite and Legendre Chaos expansions.

Module Requirements
-------------------
### Python 2.7:

- Numpy (1.13.1)

- Scipy (0.18.1)


Module Capabilities
-------------------

The `chaos_basispy` module focuses on computing the rotation using several techniques, namely:

- Gradient based method (active subspace)

- Gaussian adaptation

- Quadratic adaptation 

At this stage, the module supports classes that compute the rotation matrix that can be applied to the random variables with the respect to which the PCE is built and therefore transform them to the new random input. 

Further capabilities of computing the new coefficients of adapted PCE's by linking to application-specific forward models will be developed soon.

Demos
-----

The `demo0` and `demo1` directories contain scripts that test the performance of PolyBasis, PolyChaos and Quadrature classes that construct uni- and multi-dimensional polynomials, generate quadrature rules and compute the chaos coefficients. 
More demos on the basis adaptation techniques are coming soon...
<!--To quickly validate that the package is properly installed and bug-free, please see the /demos directory. Currently only demo1 is complete. For details on the demo please see section 3.1 in [4]. Demo2 and more are coming soon...-->

Citation
--------
      
    @misc{chaos_basispy2018,
      author = {Tsilifis, P.},
      title = {\texttt{chaos\char`_basispy}: A Polynomial Chaos basis reduction framework in python},
      howpublished = {\url{https://github.com/tsilifis/chaos_basispy}},
      year = {since 2017} 
    }

References
----------

[1] R. Tipireddy and R. Ghanem, **Basis adaptation in homogeneous chaos spaces**. *Journal of Computational Physics*, 259, pp.304-317, https://doi.org/10.1016/j.jcp.2013.12.009 (2014).

[2] P. Tsilifis and R. Ghanem, **Reduced Wiener Chaos representation of random fields via basis adaptation and projection**. *Journal of Computational Physics*, 341, pp. 102-120, https://doi.org/10.1016/j.jcp.2017.04.009 (2017).

[3] C. Thimmisetty, P. Tsilifis and R.G. Ghanem, **Polynomial Chaos basis adaptation for design optimization under uncertainty: Application to the oil well placement problem**. *Artificial Intelligence for Engineering Design, Analysis and Manufacturing*, 31(3), 265-276 https://doi.org/10.1017/S0890060417000166 (2017).

[4] P.A. Tsilifis, **Gradient-Informed Basis Adaptation for Legendre Chaos Expansions**. *ASME. J. Verif. Valid. Uncert. Quant.*, 3(1) 011005, https://doi.org/10.1115/1.4040802 (2018).

[5] P. Tsilifis and R. Ghanem, **Bayesian adaptation of chaos representations using variational inference and sampling on geodesics**. *Proc. R. Soc. A* 474 20180285, https://dx.doi.org/10.1098/rspa.2018.0285 (2018). 

[6] P. Tsilifis, X. Huan, C. Safta, K. Sargsyan, G. Lacaze, J. Oefelein, H. Najm and R. Ghanem, **Compressive sensing adaptation for polynomial chaos expansions**. *Journal of Computational Physics*, 380, pp. 29-47, https://doi.org/10.1016/j.jcp.2018.12.010 (2019).

