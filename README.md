# chaos_basispy
Polynomial Chaos Basis Reduction
=================================

Author:       Panagiotis Tsilifis

Affiliation:  Sonny Astani Department of Civil Engineering, 
              University of Southern California, Los Angeles, CA 90089, USA

email address: tsilifis@usc.edu

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
### Python :

- Numpy (1.13.1)

- Scipy (0.18.1)

*** UPDATE *** Dependence on py-orthpol has been removed and new polynomial chaos functions will be included in the current package under the module _poly_basis.


### Other :

- Dakota (only necessary for generating sparse grids, available at https://dakota.sandia.gov/) 

*** UPDATE *** New quadrature rules are available under the module _quadrature.  

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

To quickly validate that the package is properly installed and bug-free, please see the /demos directory. Currently only demo1 is complete. For details on the demo please see section 3.1 in [4]. Demo2 and more are coming soon...

Citation
--------

    @Misc{chaos_basispy2017,
      author =   {Tsilifis, P.},
      title =    {chaos_basispy: A Polynomial Chaos basis reduction framework in python},
      howpublished = {\url{https://github.com/tsilifis/chaos_basispy}},
      year = {since 2017}
    }

References
----------

[1] Tipireddy, R. and Ghanem, R., Basis adaptation in homogeneous chaos spaces. Journal of Computational Physics, 259, pp.304-317 (2014).

[2] Tsilifis, P. and Ghanem, R., Reduced Wiener Chaos representation of random fields via basis adaptation and projection. Journal of Computational Physics, 341, pp. 102-120 (2017).

[3] C. Thimmisetty, P. Tsilifis and R.G. Ghanem, Polynomial Chaos basis adaptation for design optimization under uncertainty: Application to the oil well placement problem. Artificial Intelligence for Engineering Design, Analysis and Manufacturing, 31(3), 265-276 (2017).

[4] Tsilifis, P.A., Reduced-dimensionality Legendre Chaos expansions via basis adaptation on 1d active subspaces. arXiv preprint arXiv:1611.02754 (2017).
