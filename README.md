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
to a uniform r.v through its own cdf and therefore the Legendre adapted PCE can be constructed.
Due to the above, this module considers only Hermite and Legendre Chaos expansions.

Module Capabilities
-------------------
To be written soon...
