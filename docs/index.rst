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
