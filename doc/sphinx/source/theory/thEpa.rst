.. _theoryEPA:

Electron-phonon averaging (EPA)
===============================

In Phoebe, we implemented the electron-phonon average method discuss in ref. https://doi.org/10.1016/j.mtphys.2018.07.001.
The purpose of this method is to have a fast tool for the computation of electronic transport coefficients limited by electron-phonon scattering.
The benefits of this methodology are that EPA:
1. avoids Wannierization: while accurate, finding Wannier functions is tricky and is not yet fully automatizable; and
2. gives priority to speed rather than accuracy, allowing a bigger computational throughput.

We only consider transport coefficients within the relaxation time approximation.
In such case, the electronic lifetimes due to the electron-phonon scattering are:

.. math::

   \tau_{n\boldsymbol{k}}^{-1}(\mu,T)
   =
   \frac{V}{(2\pi)^2 \hbar} \sum_{m\nu}
   \int_{BZ} d\boldsymbol{q}
   |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
   \bigg[ \big(n(\omega_{\boldsymbol{q}\nu},T) + f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} + \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})  \\
   + \big(n(\omega_{\boldsymbol{q}\nu},T) + 1 - f(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m},\mu,T)\big) \delta(\epsilon_{\boldsymbol{k}n} - \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}+\boldsymbol{q}m}) \bigg]


While accurate, this expression is still too expensive for some purposes.


The main idea of EPA is to replace integrals over the Brillouin zone with some averaged quantities.
In particular, rather than working with the full phonon dispersion, we average the frquencies of each phonon mode:

.. math::

   \omega_{\boldsymbol{q}\nu}
   \to
   \bar{\omega}_{\nu}

So that we have :math:`3 N_{atoms}` modes per crystal.
Similarly, the make an energy binning on the electron-phonon coupling strength, as

.. math::

   |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2
   \to
   g^2_{\nu} (\epsilon_{\boldsymbol{k}n}, \epsilon_{\boldsymbol{k}+\boldsymbol{q}m})

There can be various ways of doing this averaging procedure on :math:`g`.
In Phoebe, we implemented a gaussian quadrature procedure.
In this approach, the electron-phonon is approximated as a weighted sum of the electron-phonon coupling computed from an ab-initio code:

.. math::
   g^2_{\nu} (\epsilon_1,\epsilon_2)
   =
   \frac{1}{W}
   \sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}} |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2

where

.. math::

   W = \sum_{mn\boldsymbol{k}\boldsymbol{q}} w_{mn\boldsymbol{k}\boldsymbol{q}}

and the weights :math:`w` are found by minimizing

.. math::
   \sum_{mn} \int d\boldsymbol{q} d\boldsymbol{k} ( g^2_{\nu} (\epsilon_1,\epsilon_2) - |g_{mn\nu}(\boldsymbol{k},\boldsymbol{q})|^2 )^2
   \exp\bigg( -\frac{(\epsilon_{\boldsymbol{k}n}-\epsilon_1)^2+(\epsilon_{\boldsymbol{k}+\boldsymbol{q}m}-\epsilon_2)^2}{2\sigma^2_{gauss}} \bigg)

After the electron-phonon coupling is approximated in this way, the electron lifetimes can be integrated as:

.. math::
   \tau^{-1}(\epsilon,\mu,T)
   =
   \frac{2\pi V}{g_s \hbar} \sum_{\nu}
   g^2_{\nu}(\epsilon,\epsilon+\bar{\omega}_{\nu})
   \big(n(\bar{\omega}_{\nu},T) + f(\epsilon + \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon + \bar{\omega}_{\nu})  +  \\
   g^2_{\nu}(\epsilon,\epsilon-\bar{\omega}_{\nu})
   \big(n(\bar{\omega}_{\nu},T) + f(\epsilon - \bar{\omega}_{\nu},\mu,T)\big) \rho(\epsilon - \bar{\omega}_{\nu})

where

.. math::
   v^2_{\alpha\beta} (\epsilon) \rho (\epsilon)
   =
   \sum_n \int_{BZ} d\boldsymbol{k} v_{n\boldsymbol{k}\alpha} v_{n\boldsymbol{k}\beta} \delta(\epsilon-\epsilon_{\boldsymbol{k}n})

Finally, transport coefficients can be computed from the following tensor:

.. math::
   K_{\alpha\beta}^{(p)}
   =
   \frac{g_s e^{2-p}}{(2\pi)^{3} (k_B T)^{p+1}} \int d\epsilon v^2_{\alpha\beta}(\epsilon) \rho(\epsilon) \tau(\epsilon,\mu,T) (\epsilon-\mu)^p f(\epsilon,\mu,T) [1-f(\epsilon,\mu,T)]

For example, the electrical conductivity is:

.. math::
   \sigma_{\alpha\beta}(\mu,T) = K_{\alpha\beta}^{(0)}

Note that, in order to compute the velocity term, we still need a form of interpolation of the electronic band structure.
In order to keep the spirit of avoiding the Wannierization procedure, we use the @ref thFourierElectrons.

Furthermore, we note that the density of states is integrated using the tetrahedron method.

In terms of computational parameters (besides temperature and doping concentration) the EPA requires tuning a few parameters, namely

1. the size of the coarse k/q grid used in an ab-initio code to compute the coupling;

2. the energy bins used to approximate the electron-phonon coupling (which are set when converting data from Quantum ESPRESSO to Phoebe),

3. the energy bins used to integrate lifetimes (typically, one should use a better energy sampling than the one for the coupling)

4. the grid of wavevectors, used to average velocities and density of states.


Electronic Fourier interpolation
----------------------------------

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
