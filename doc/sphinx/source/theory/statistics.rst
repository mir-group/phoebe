Basic particle properties
===============================

Particle statistics
-------------------

**For electrons:** The Fermi--Dirac distribution for a Bloch state (:math:`\boldsymbol{k}`, :math:`b`) (:math:`\boldsymbol{k}`, :math:`b` are Bloch indices) represents the occupation of electron states, and is: 

.. math::
   f_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}-\mu}{k_BT})}+1}

**For phonons:** The Bose-Einstein distribution, which represents phonon occupation, is: 

.. math::
   n_{\boldsymbol{k},b} = \frac{1}{e^{(\frac{\epsilon_{\boldsymbol{k},b}}{k_BT})}-1}

The derivatives of these functions show up in many transport equations. A computationally stable way to evaluate derivatives of these occupation functions is:

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

The specific heat at constant volume for phonons (or electrons, substituting in the appropriate occupation factors) is evaluated as:

.. math::
   C_v = \frac{1}{V N_k} \sum_{\boldsymbol{k},b} \epsilon_{\boldsymbol{k},b} \frac{\partial n_{\boldsymbol{k},b}}{\partial T}


where :math:`V` is the volume of the primitive crystal unit cell, :math:`N_k` is the number of wavevectors used, :math:`\boldsymbol{k}` is a wavevector index, :math:`b` is a branch/band index, :math:`\epsilon` is the phonon energy (for electrons, this is replaced by :math:`\epsilon-\mu`), :math:`T` the temperature and :math:`n` is the Bose--Einstein (Fermi--Dirac) distribution function.



Density of States
-----------------

The density of states is defined as the number of states that are available to a particle at a certain energy, :math:`E`:

.. math::
   DOS(E) = \frac{1}{(2\pi)^d V} \sum_b \int_{BZ} \delta(E-\epsilon_{\boldsymbol{k}b}) d\boldsymbol{k} 

where :math:`d` is the dimensionality, :math:`V` is the crystal unit cell volume, :math:`b` is a Bloch index over bands, :math:`\boldsymbol{k}` is a Bloch index over wavevectors, :math:`\epsilon` is the particle energy and the integration is carried over the Brillouin zone.

The definition holds for both phonon and electrons.

The integral is sampled with a uniform grid of points in the Brillouin zone.
In the DoS apps, the integral of the Dirac delta is implemented using the tetrahedron method.


Dirac delta approximations
--------------------------

We offer two possible methods for the approximation of the Dirac delta functions used in Phoebe's transport calculations. 

.. raw:: html

  <h4>Gaussian Approximation</h4>

The delta function for the energy conservation can be replaced by a Gaussian: 

.. math::
   \delta(\hbar \omega)=\frac{1} {\sqrt{\pi}  \sigma} \exp{\left[-(\hbar \omega/ \sigma )^2 \right]} \;,

where the width, :math:`\sigma` is a constant decided by user input.
It is important to note that when the delta function is substituted with a Gaussian, the detailed balance condition is only valid under approximation.
The definition used above guarantees that the scattering matrix is symmetric and non-negative.

.. raw:: html

  <h4>Adaptive Gaussian</h4>

Another method is the adaptive Gaussian smearing scheme (see https://link.aps.org/doi/10.1103/PhysRevB.75.195121).
We use this method to approximate a Dirac delta function of the form:

.. math::
   \delta(\hbar (\omega_1+\omega_2-\omega_3)=\frac{1} {\sqrt{\pi}  \sigma} \exp{(-(\hbar (\omega_1+\omega_2-\omega_3)/ \sigma )^2)}.

In this method, :math:`\sigma` is now dependent on the energies, and is not a user input.
Specifically, we build it as:

.. math::
   \sigma = \frac{1}{\sqrt{12}} \sqrt{ \sum_{\beta} \left(\sum_{\alpha} (v_2-v_3\right) \frac{M_{\alpha \beta}}{N_{\beta}}  )^2 }

where :math:`M` is a matrix comprised of the primitive cell lattice vectors (each column is a lattice vector), :math:`v_2` and :math:`v_3` are phonon group velocities, and :math:`N_{\beta}` is the number of wavevectors sampled along direction :math:`\beta`.

Note that the adaptive scheme may be critical in the case where the velocity sum to zero: in that case, we skip the scattering event, unless we have an exact energy conservation taking place.


Dynamical matrix
-----------------

A density-functional code, detailed elsewhere, can compute the following force-constants matrix:

.. math::
   M(ls\alpha | l's'\alpha') = \frac{\partial^2 \mathcal{E}}{\partial u_{ls\alpha} \partial u_{l's'\alpha'}}

where :math:`M` is a matrix of second order derivative of the total crystal energy :math:`\mathcal{E}` with respect to the ionic displacement :math:`u_{ls\alpha}`, where :math:`l` labels a unit cell in a supercell, :math:`s` is an index over the ionic basis, and :math:`\alpha` denotes the direction in which the displacement is made.
This matrix can either be computed with density functional perturbation theory or with a frozen-phonon approach.
Due to the periodicity of the crystal, one can set :math:`l=0`.

The dynamical matrix is the Fourier transform of this matrix.
Excluding polar corrections, the dynamical matrix is:

.. math::
   D(s\alpha | s'\alpha')(\boldsymbol{q}) = \sum_{l'} M(0s\alpha | l's'\alpha') e^{i \boldsymbol{q} \cdot \boldsymbol{R}_{l'}}

Note that the Bravais lattice vectors are defined as the Bravais lattice vectors belonging to the Wigner-Seitz zone (not the Brillouin zone!) of a supercell, whose size is :math:`N_{qx}\times N_{qy}\times N_{qz}` that of the primitive unit cell and this is the size of the q-point mesh used to compute the phonons in the DFT code.
   
If ions carry a charge, one must not forget to add an additional term to D:

.. math::
   D(s\alpha | s'\alpha')(\boldsymbol{q}) += \frac{4\pi}{\Omega} e^2 \frac{ (\boldsymbol{q} \cdot Z^*_s)_{\alpha} (\boldsymbol{q} \cdot Z^*_{s'})_{\alpha'} } { (\boldsymbol{q} \cdot \epsilon^{\infty} \cdot \boldsymbol{q}) }

where :math:`Z_{s,\alpha,\beta}` is the Born charge tensor of atom :math:`s` and :math:`\epsilon^{\infty}` is the static dielectric constant.

The phonon energy and phonon eigenvectors are defined from the diagonalization problem as:

.. math::
   D(s\alpha | s'\alpha')(\boldsymbol{q}) z_{s'\alpha'j}(\boldsymbol{q}) = \omega_{j}^2(\boldsymbol{q}) z_{s\alpha j}(\boldsymbol{q})



