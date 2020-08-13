@page Theory Theory
@section MAT2 Dynamical matrix

A density-functional code, detailed elsewhere, can compute the following force-constants matrix:
\begin{equation}
M(ls\alpha | l's'\alpha') = \frac{\partial^2 \mathcal{E}}{\partial u_{ls\alpha} \partial u_{l's'\alpha'}}
\end{equation}
where \f$M\f$ is a matrix of second order derivative of the total crystal energy \f$\mathcal{E}\f$ with respect to the ionic displacement \f$u_{ls\alpha}\f$, where \f$l\f$ labels a unit cell in a supercell, \f$s\f$ is an index over the ionic basis, and \f$\alpha\f$ denotes the direction in which the displacement is made.
This matrix can either be computed with density functional perturbation theory or with a frozen-phonon approach.
Due to the periodicity of the crystal, one can set \f$l=0\f$.

The dynamical matrix is the Fourier transform of this matrix.
Excluding polar corrections, the dynamical matrix is:
\begin{equation}
D(s\alpha | s'\alpha')(q) = \sum_{l'} M(0s\alpha | l's'\alpha') e^{i q \cdot R_{l'}}
\end{equation}
If ions carry a charge, one must not forget to add an additional term to D:
\begin{equation}
D(s\alpha | s'\alpha')(q) += \frac{4\pi}{\Omega} e^2 \frac{ (q \cdot Z^*_s)_{\alpha} (q \cdot Z^*_{s'})_{\alpha'} } { (q \cdot \epsilon^{\infty} \cdot q) }
\end{equation}
where \f$Z_{s,\alpha,\beta}\f$ is the Born charge tensor of atom \f$s\f$ and \f$\epsilon^{\infty}\f$ is the static dielectric constant.

The phonon energy and phonon eigenvectors are defined from the diagonalization problem as:

\begin{equation}
D(s\alpha | s'\alpha')(q) z_{s'\alpha'j}(q) = \omega_{j}^2(q) z_{s\alpha j}(q)
\end{equation}

