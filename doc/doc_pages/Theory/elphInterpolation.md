@page Theory Theory
@section ELPHC Electron-phonon coupling

This is a summary of the strategy for interpolation adopted by Phoebe, as detailed also at this [reference] (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.165108).

The coupling to be used for the calculation of scattering rates is:


\begin{equation}
g^{SE}_{mn,\nu} (k,q) = \bigg( \frac{1}{2 m \omega_{q\nu}} \bigg)^{1/2} g_{mn,\nu} (k,q)
\end{equation}

What we need to interpolate instead is:
\begin{equation}
g_{mn,\nu} (k,q) = \big<mk+q \big| \partial_{q\nu}V \big| nk \big>
\end{equation}

First recall the relation between the wavefunction in Wannier and Bloch representation.

\begin{equation}
\big|mR_e\big> = \sum_{nk} e^{-ikR_e} U_{nm,k} \big|nk\big>
\end{equation}

\begin{equation}
\big|nk\big> = \frac{1}{N_e} \sum_{mR_e} e^{ikR_e} U_{mn,k}^\dagger \big|mR_e\big>
\end{equation}
where \f$N_e\f$ is the number of supercells.



Let \f$e_{q\kappa}^{\nu}\f$ be the phonon eigenvector we get from diagonalizing the dynamical matrix.
We define \f$u_{q\kappa}^{\nu} = (\frac{m_0}{m_{\kappa}})^{1/2} e_{q\kappa}^{\nu}\f$, where \f$m_0\f$ is the electron mass, and \f$m_{\kappa}\f$ is the ionic mass.

To transform the potential from the reciprocal to the real space representation, we have:
\begin{equation}
\partial_{\kappa R_p} V(r)
=
\frac{1}{N_p}
\sum_{q\nu} e^{-iqR_p} [u_{q\kappa}^{\nu}]^{-1} \partial_{q\nu} V(r)
\end{equation}


So, we first transform to Wannier space by:
\begin{equation}
g(R_e,R_p)
=
\frac{1}{N_p}
\sum_{kq} e^{-ikR_e-iqR_p} U_{k+q}^\dagger g(k,q) U_k u_q^{-1}
\end{equation}
Currently, the code phoebe is reading the matrix \f$g(R_e,R_p)\f$ from the code EPW.

Then, we interpolate to Bloch space

\begin{equation}
g(k,q)
=
\frac{1}{N_e}
\sum_{R_e R_p} e^{ikR_e+iqR_p} U_{k+q} g(R_e,R_p) U_k^\dagger u_q
\end{equation}


Details:
* the mesh of k and q points must be the same (or at least commensurate, so that we can map k+q into the same k grid).

* grids must be unshifted, so that k+q falls into the same grid.
Also, if k+q is outside the 1st BZ, we can use the parallelt transport gauge and set \f$U_{k+q+G} = U_{k+q}\f$.

* q computed by an ab-initio code is in the irreducible wedge, while k is on a complete uniform grid.

