Electronic Fourier interpolation
================================

We implemented the Fourier interpolation of the band structure (also used in Boltztrap v1).  
The method is well described [in this reference] (https://link.aps.org/doi/10.1103/PhysRevB.38.2721).
Note that the algorithm described below works for a single band (so it must be repeated for every distinct band we want to interpolate).

Let's suppose to have :math:`N` data points, i.e. an energy :math:`\epsilon(\boldsymbol{k}_i)`, specified over a (coarse) mesh of kpoints :math:`\boldsymbol{k}_i`, with :math:`i=0,\dots,N-1`.
We want to interpolate these points, so we can obtain an energy :math:`\tilde{\epsilon}(\boldsymbol{k})` for an arbitrary k-point.

We define the interpolating function as:

.. math::
   \tilde{\epsilon}(\boldsymbol{k}) = \sum_{m=0}^{M-1} c_m S_m(\boldsymbol{k}) \;,

where :math:`c_m` are expansion coefficients (to be found) and

.. math::
   S_m(\boldsymbol{k}) = \frac{1}{n} \sum_{\Lambda} e^{i\boldsymbol{k} \Lambda \boldsymbol{R}_m} \;,

is a star function, where :math:`\Lambda` is a point-group symmetry operation of the crystal, :math:`n` is the number of symmetry operations, and :math:`\boldsymbol{R}_m` is a lattice vector.

The choice of :math:`\boldsymbol{R}_m` is a free parameter of the interpolation algorithm, and the user can fix it by providing a cutoff, identifying all :math:`\boldsymbol{R}_m` such that :math:`|\boldsymbol{R}_m | < R_{\text{cut}}`.
We label vectors such as :math:`m=0,\dots,M-1`, and :math:`m=0` identifies the null vector.
Note that one must provide more lattice vectors than points available in the system.

To find the expansion coefficients, we minimize a Lagrangian :math:`\mathcal{L}` under the constraint that the function interpolates the data points.
In particular, we want to minimize:

.. math::
   \mathcal{L} = \frac{1}{2} \sum_m c_m \rho_m + \sum_i \lambda_i (\epsilon(\boldsymbol{k}_i)-\tilde{\epsilon}(\boldsymbol{k}_i)) \;,

where :math:`\lambda_i` is a set of Lagrange multipliers, and the roughness function :math:`\rho_m` is defined as:

.. math::
   \rho_m = \bigg(1-A\frac{R_m}{R_{min}}\bigg)^2 + B\bigg(\frac{R_m}{R_{min}}\bigg)^6  \;,

where we fix :math:`A=B=3/4`, and :math:`R_{min}` is the norm of the smallest non-zero lattice vector.

After solving the Lagrange problem, one can compute the Lagrange multipliers from the following linear algebra problem.
Choose a particular reference point, here we use i=0.
Construct the matrix :math:`H` (of size N-1) as

.. math::
   H_{ij} = \sum_{m=1}^{M-1} \frac{ (S_m(\boldsymbol{k}_i)-S_m(\boldsymbol{k}_0)) (S^*_m(\boldsymbol{k}_j) - S^*_m(\boldsymbol{k}_0)) }{\rho_m} \;,

and solve the linear algebra problem:

.. math::
   \sum_{j=1}^{N-1} H_{ij} \lambda_j = \epsilon_{k_i} - \epsilon_{k_0} \;.

After having obtained the Lagrange multipliers, the expansion coefficients are found as:

.. math::
   c_m = \rho_m^{-1} \sum_{i=1}^{N-1} \lambda_i ( S^*_m(\boldsymbol{k}_i) - S^*_m(\boldsymbol{k}_0) ) \;,

and for the zero lattice vector:

.. math::
   c_0 = \epsilon(\boldsymbol{k}_0) - \sum_{m=1}^{M-1} c_m S_m(k_0) \;.

One can compute the expansion coefficients once and store them in memory.
The star function :math:`S` must be recomputed at every evaluation of energy.
Additionally, the velocity is easily computed as:

.. math::
   \tilde{v}(\boldsymbol{k}) = \sum_{m=0}^{M-1} c_m \bigg( i \frac{1}{n} \sum_{\Lambda}  \Lambda \boldsymbol{R}_m e^{i\boldsymbol{k} \Lambda \boldsymbol{R}_m} \bigg) \;.
