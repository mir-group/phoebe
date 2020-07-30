@page Theory
@section FINTERP Fourier interpolation of the electronic band structure

We implemented the Fourier interpolation of the band structure (also used in Boltztrap v1).  
The method is well described in Ref. [Pickett, Krakauer and Allen; PRB 38, 2721 (1988)].
Note that the algorithm works for a single band (so it must be repeated for every distinct band we want to interpolate).

Let's suppose to have \f$N\f$ data points, i.e. an energy \f$\epsilon(\boldsymbol{k}_i)\f$, specified over a (coarse) mesh of kpoints \f$\boldsymbol{k}_i\f$, with \f$i=0,\dots,N-1\f$.
We want to interpolate these points, so we can obtain an energy \f$\tilde{\epsilon}(\boldsymbol{k})\f$ for an arbitrary k-point.

We define the interpolating function as:
\begin{equation}
\tilde{\epsilon}(\boldsymbol{k}) = \sum_{m=0}^{M-1} c_m S_m(\boldsymbol{k}) \;,
\end{equation}
where $c_m$ are expansion coefficients (to be found) and
\begin{equation}
S_m(\boldsymbol{k}) = \frac{1}{n} \sum_{\Lambda} e^{i\boldsymbol{k} \Lambda \boldsymbol{R}_m} \;,
\end{equation}
is a star function, where $\Lambda$ is a point-group symmetry operation of the crystal, \f$n\f$ is the number of symmetry operations, and \f$\boldsymbol{R}_m\f$ is a lattice vector.

The choice of \f$\boldsymbol{R}_m\f$ is a free parameter of the interpolation algorithm, and the user can fix it by providing a cutoff, identifying all \f$\boldsymbol{R}_m\f$ such that \f$|\boldsymbol{R}_m | < R_{\text{cut}}\f$.
We label vectors such as \f$m=0,\dots,M-1\f$, and \f$m=0\f$ identifies the null vector.
Note that one must provide more lattice vectors than points available in the system.

To find the expansion coefficients, we minimize a Lagrangian \f$\mathcal{L}\f$ under the constraint that the function interpolates the data points.
In particular, we want to minimize:
\begin{equation}
\mathcal{L} = \frac{1}{2} \sum_m c_m \rho_m + \sum_i \lambda_i (\epsilon(\boldsymbol{k}_i)-\tilde{\epsilon}(\boldsymbol{k}_i)) \;,
\end{equation}
where $\lambda_i$ is a set of Lagrange multipliers, and the roughness function $\rho_m$ is defined as:
\begin{equation}
\rho_m = \bigg(1-A\frac{R_m}{R_{min}}\bigg)^2 + B\bigg(\frac{R_m}{R_{min}}\bigg)^6  \;,
\end{equation}
where we fix \f$A=B=3/4\f$, and \f$R_{min}\f$ is the norm of the smallest non-zero lattice vector.

After solving the Lagrange problem, one can compute the Lagrange multipliers from the following linear algebra problem.
Choose a particular reference point, here we use i=0.
Construct the matrix \f$H\f$ (of size N-1) as
\begin{equation}
H_{ij} = \sum_{m=1}^{M-1} \frac{ (S_m(\boldsymbol{k}_i)-S_m(\boldsymbol{k}_0)) (S^*_m(\boldsymbol{k}_j) - S^*_m(\boldsymbol{k}_0)) }{\rho_m} \;,
\end{equation}
and solve the linear algebra problem:
\begin{equation}
\sum_{j=1}^{N-1} H_{ij} \lambda_j = \epsilon_{k_i} - \epsilon_{k_0} \;.
\end{equation}

After having obtained the Lagrange multipliers, the expansion coefficients are found as:
\begin{equation}
c_m = \rho_m^{-1} \sum_{i=1}^{N-1} \lambda_i ( S^*_m(\boldsymbol{k}_i) - S^*_m(\boldsymbol{k}_0) ) \;,
\end{equation}
and for the zero lattice vector:
\begin{equation}
c_0 = \epsilon(\boldsymbol{k}_0) - \sum_{m=1}^{M-1} c_m S_m(k_0) \;.
\end{equation}

One can compute the expansion coefficients once and store them in memory.
The star function \f$S\f$ must be recomputed at every evaluation of energy.
Additionally, the velocity is easily computed as:
\begin{equation}
\tilde{\epsilon}(\boldsymbol{k}) = \sum_{m=0}^{M-1} c_m \bigg( i \frac{1}{n} \sum_{\Lambda}  \Lambda \boldsymbol{R}_m e^{i\boldsymbol{k} \Lambda \boldsymbol{R}_m} \bigg) \;.
\end{equation}


