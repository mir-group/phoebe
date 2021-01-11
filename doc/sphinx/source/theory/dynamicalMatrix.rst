Dynamical matrix
================

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


