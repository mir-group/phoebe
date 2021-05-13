Independent particle properties
===============================

Particle statistics
-------------------

Equations available in class Particle.
The Fermi--Dirac distribution (for electrons) for a Bloch state (:math:`\boldsymbol{k}`, :math:`b`) (:math:`\boldsymbol{k}`, :math:`b` are Bloch indices) is

.. math::
   f_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}-\mu}{k_BT})}+1}

The Bose--Einstein distribution (for phonons) is:

.. math::
   n_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}}{k_BT})}-1}

A computationally stable way to evaluate derivatives is:

.. math::
   \frac{\partial n_{\boldsymbol{k},b}}{\partial T} = \frac{\epsilon_{\boldsymbol{k},b}}{4k_BT^2} \sinh( \frac{\epsilon_{\boldsymbol{k},b}}{2T} )

.. math::
   \frac{\partial f_{\boldsymbol{k},b}}{\partial T} = \frac{\epsilon_{\boldsymbol{k},b}-\mu}{4k_BT^2} \cosh( \frac{\epsilon_{\boldsymbol{k},b}-\mu}{2T} )

.. math::
   \frac{\partial n_{\boldsymbol{k},b}}{\partial \epsilon} = - \frac{1}{4k_BT} \sinh( \frac{\epsilon_{\boldsymbol{k},b}}{2T} )

.. math::
   \frac{\partial f_{\boldsymbol{k},b}}{\partial \epsilon} = - \frac{1}{4k_BT} \cosh( \frac{\epsilon_{\boldsymbol{k},b}-\mu}{2T} )


Specific heat
-------------

The specific heat at constant volume for phonons (electrons) is evaluated as:

.. math::
   C_v = \frac{1}{V N_k} \sum_{\boldsymbol{k},b} \epsilon_{\boldsymbol{k},b} \frac{\partial n_{\boldsymbol{k},b}}{\partial T}


where :math:`V` is the volume of the primitive crystal unit cell, :math:`N_k` is the number of wavevectors used, :math:`\boldsymbol{k}` is a wavevector index, :math:`b` is a branch/band index, :math:`\epsilon` is the phonon(electron, w.r.t. the chemical potential) energy, :math:`T` the temperature and :math:`n` is the Bose--Einstein (Fermi--Dirac) distribution function.





Density of States
-----------------

The density of states is defined as the number of states that are available to a particle at a certain energy :math:`E`:

.. math::
   DOS(E) = \frac{1}{(2\pi)^d V} \sum_b \int_{BZ} d\boldsymbol{k} \delta(E-\epsilon_{\boldsymbol{k}b})

where :math:`d` is the dimensionality, :math:`V` is the crystal unit cell volume, :math:`b` is a Bloch index over bands, :math:`\boldsymbol{k}` is a Bloch index over wavevectors, :math:`\epsilon` is the particle energy and the integration is carried over the Brillouin zone.

The definition holds for both phonon and electrons.

The integral is sampled with a uniform grid of points in the Brillouin zone.
In the DoSApps, the integral of the Dirac-delta is implemented using the tetrahedron method.
