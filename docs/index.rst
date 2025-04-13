.. module:: cosmology.compat.camb

Cosmology API compatibility layer for CAMB
==========================================

The :mod:`cosmology.compat.camb` package provides a compatibility layer between
CAMB and the `Cosmology API <https://cosmology.readthedocs.io>`_.


Usage
-----

To create a Cosmology API-compliant ``cosmo`` object, wrap CAMB's ``pars`` and
``results`` in ``cosmology.compat.camb.Cosmology``:

.. code-block::

   import camb
   import cosmology.compat.camb

   pars = camb.set_params(H0=70.0)
   results = camb.get_background(pars)

   cosmo = cosmology.compat.camb.Cosmology(pars, results)


Reference
---------

.. currentmodule:: cosmology.compat.camb

.. autodata:: K_LINEAR

.. autoproperty:: Cosmology.h
.. autoproperty:: Cosmology.H0
.. automethod:: Cosmology.H
.. automethod:: Cosmology.H_over_H0

.. autoproperty:: Cosmology.Omega_m0
.. automethod:: Cosmology.Omega_m
.. autoproperty:: Cosmology.Omega_de0
.. automethod:: Cosmology.Omega_de
.. autoproperty:: Cosmology.Omega_k0
.. automethod:: Cosmology.Omega_k

.. autoproperty:: Cosmology.critical_density0

.. autoproperty:: Cosmology.hubble_distance
.. automethod:: Cosmology.comoving_distance
.. automethod:: Cosmology.transverse_comoving_distance
.. automethod:: Cosmology.angular_diameter_distance

.. automethod:: Cosmology.inv_comoving_distance

.. automethod:: Cosmology.growth_factor
