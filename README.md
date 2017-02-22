# chaos_basispy
Polynomial Chaos Basis Adaptation
=================================

Author:       Panagiotis Tsilifis

Affiliation:  Sonny Astani Department of Civil Engineering, 
              University of Southern California, Los Angeles, CA 90089, USA


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

- Numpy

- Scipy

- Py-orthpol (available at https://github.com/PredictiveScienceLab/py-orthpol)

### Other :

- Dakota (only necessary for generating sparse grids, available at https://dakota.sandia.gov/) 

Module Capabilities
-------------------

The `chaos_basispy` module focuses on computing the rotation using several techniques, namely:

- Active subspace method

- Gaussian adaptation

- Quadratic adaptation 

More detailed description to be written soon...

Demos
-----

To quickly validate that the package is properly installed and bug-free, please see the /demos directory. Currently only demo1 is complete. For details on the demo please see section 3.1 in [4]. Demo2 and more are coming soon...

Citation
--------

    @Misc{chaos_basispy2017,
      author =   {Tsilifis, P.},
      title =    {chaos_basispy: A Polynomial Chaos basis adaptation framework in python},
      howpublished = {\url{https://github.com/tsilifis/chaos_basispy}},
      year = {since 2017}
    }

References
----------

[1] C. Thimmisetty, P. Tsilifis and R.G. Ghanem, Polynomial Chaos basis adaptation for design optimization under uncertainty: Application to the oil well placement problem. To appear in Artificial Intelligence for Engineering Design, Analysis and Manufacturing, 2016.

[2] Tipireddy, R. and Ghanem, R., Basis adaptation in homogeneous chaos spaces. Journal of Computational Physics, 259, pp.304-317, 2014.

[3] Tsilifis, P. and Ghanem, R., Reduced Wiener Chaos representation of random fields via basis adaptation and projection. To appear in Journal of Computational Physics (arXiv:1603.04803), 2016.

[4] Tsilifis, P.A., Reduced-dimensionality Legendre Chaos expansions via basis adaptation on 1d active subspaces. arXiv preprint arXiv:1611.02754, 2016.
