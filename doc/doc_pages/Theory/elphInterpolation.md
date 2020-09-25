@page Theory Theory
@section ELPHC Electron-phonon coupling, Wannier interpolation

This is a summary of the strategy for interpolation adopted by Phoebe, as detailed also at this [reference] (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.165108).

The coupling to be used for the calculation of scattering rates is:


\begin{equation}
g^{SE}_{b'b,\nu} (\boldsymbol{k},\boldsymbol{q}) = \bigg( \frac{1}{2 m \omega_{\boldsymbol{q}\nu}} \bigg)^{1/2} g_{b'b,\nu} (\boldsymbol{k},\boldsymbol{q})
\end{equation}

What we need to interpolate instead is:
\begin{equation}
g_{b'b,\nu} (\boldsymbol{k},\boldsymbol{q}) = \big<b'\boldsymbol{k}+\boldsymbol{q} \big| \partial_{\boldsymbol{q}\nu}V \big| b\boldsymbol{k} \big>
\end{equation}

First recall the relation between the wavefunction in Wannier and Bloch representation.

\begin{equation}
\big|m\boldsymbol{R}_e\big> = \sum_{b\boldsymbol{k}} e^{-i\boldsymbol{k}\cdot\boldsymbol{R}_e} U_{mb,\boldsymbol{k}} \big|b\boldsymbol{k}\big>
\end{equation}

\begin{equation}
\big|b\boldsymbol{k}\big> = \frac{1}{N_e} \sum_{m\boldsymbol{R}_e} e^{i\boldsymbol{k}\cdot\boldsymbol{R}_e} U_{bm,\boldsymbol{k}}^\dagger \big|m\boldsymbol{R}_e\big>
\end{equation}
where \f$N_e\f$ is the number of supercells.



Let \f$e_{\boldsymbol{q}\kappa}^{\nu}\f$ be the phonon eigenvector we get from diagonalizing the dynamical matrix.
We define \f$u_{\boldsymbol{q}\kappa}^{\nu} = (\frac{m_0}{m_{\kappa}})^{1/2} e_{\boldsymbol{q}\kappa}^{\nu}\f$, where \f$m_0\f$ is the electron mass, and \f$m_{\kappa}\f$ is the ionic mass.

To transform the potential from the reciprocal to the real space representation, we have:
\begin{equation}
\partial_{\kappa \boldsymbol{R}_p} V(\boldsymbol{r})
=
\frac{1}{N_p}
\sum_{\boldsymbol{q}\nu} e^{-i\boldsymbol{q}\cdot\boldsymbol{R}_p} [u_{\boldsymbol{q}\kappa}^{\nu}]^{-1} \partial_{\boldsymbol{q}\nu} V(\boldsymbol{r})
\end{equation}


So, we first transform to Wannier space by:
\begin{equation}
g(\boldsymbol{R}_e,\boldsymbol{R}_p)
=
\frac{1}{N_p}
\sum_{\boldsymbol{k}\boldsymbol{q}} e^{-i\boldsymbol{k}\cdot\boldsymbol{R}_e-i\boldsymbol{q}\cdot\boldsymbol{R}_p} U_{\boldsymbol{k}+\boldsymbol{q}}^\dagger g(\boldsymbol{k},\boldsymbol{q}) U_{\boldsymbol{k}} u_{\boldsymbol{q}}^{-1}
\end{equation}
Currently, the code phoebe is reading the matrix \f$g(\boldsymbol{R}_e,\boldsymbol{R}_p)\f$ from the code EPW.

Then, we interpolate to Bloch space

\begin{equation}
g(\boldsymbol{k},\boldsymbol{q})
=
\frac{1}{N_e}
\sum_{\boldsymbol{R}_e \boldsymbol{R}_p} e^{i\boldsymbol{k}\cdot\boldsymbol{R}_e+i\boldsymbol{q}\cdot\boldsymbol{R}_p} U_{\boldsymbol{k}+\boldsymbol{q}} g(\boldsymbol{R}_e,\boldsymbol{R}_p) U_{\boldsymbol{k}}^\dagger u_{\boldsymbol{q}}
\end{equation}


Details:
* the mesh of \f$\boldsymbol{k}\f$ and \f$\boldsymbol{q}\f$ points must be the same (or at least commensurate, so that we can map \f$\boldsymbol{k}+\boldsymbol{q}\f$ into the same \f$\boldsymbol{k}\f$ grid).






@subsection COMPELPH Computational details, Wannier90

First, the mesh of \f$\boldsymbol{k}\f$ and \f$\boldsymbol{q}\f$ points must be commensurate.
In fact, the matrices \f$U\f$ are computed on a grid of \f$\boldsymbol{k}\f$ points.
It is thus necessary that \f$\boldsymbol{k}+\boldsymbol{q}\f$ falls on the same grid of points.


To be precise, note that we use two different sets of \f$U\f$ matrices.
When interpolating from the Wannier to the Bloch space, we use the U matrices computed at an arbitrary grid of points obtained by diagonalizing the Hamiltonian matrix in the Wannier space.
When first transforming from the Bloch space (the coupling with an ab-initio code) to the Wannier space, we use the \f$U\f$ matrices computed by Wannier90.
This is necessary, because the coupling \f$g\f$ is a complex function, and we must rotate it using the same wavefunction gauge used to build the maximally localized wannier functions.
Moreover, Wannier90 may use the disentanglement procedure.
In that case, the Bloch to Wannier transformation is:
\begin{equation}
\big|m\boldsymbol{R}_e\big> = \sum_{\boldsymbol{k} b b'} e^{-i\boldsymbol{k}\cdot\boldsymbol{R}_e} U^{rot}_{mb',\boldsymbol{k}} U^{dis}_{b'b,\boldsymbol{k}} \big|b\boldsymbol{k}\big>
\end{equation}
where the number of disentangled bands \f$b'\f$ is smaller than the number of entangled bands \f$b\f$.
Therefore, we rotate the electron-phonon coupling from the Bloch to Wannier space using the entangled number of bands.
Wannier90 prints the two different \f$U\f$ matrices, and one can just multiply them to get the transformation matrix.

As a further minor detail, remember that some bands (like deep core bands) may be excluded from the Wannierization procedure (through the keyword exclude-indeces), so that there may be an offset in the band index of U and g.



\subsection COMPELPHQE Computational details, gauge fixing in Quantum ESPRESSO 

The interpolation procedure described above implicitely assumes that the wavefunction \f$\big|b\boldsymbol{k}\big>\f$ has a fixed gauge.
This is because the rotation \f$U\f$ must act on the same wavefunction used to compute the coupling \f$g\f$.

If no gauge is imposed on the Bloch wavefunction, \f$g\f$ is a random function in the Brillouin zone.
In details: consider the Bloch Hamiltonian \f$H_{\boldsymbol{k}}\f$.
A DFT code using periodic boundary conditions will diagonalize it to find \f$ H_{\boldsymbol{k}} \psi_{\boldsymbol{k}} = \epsilon_k \psi_{\boldsymbol{k}}\f$.
However, we can adopt the transformation \f$\psi_{\boldsymbol{k}} \to e^{i \theta_{\boldsymbol{k}}} \psi_{\boldsymbol{k}}\f$ and still have \f$e^{i \theta_{\boldsymbol{k}}} \psi_{\boldsymbol{k}}\f$ be an eigenvector.
Specifically, the numerical diagonalization procedure doesn't by itself fix a gauge: in fact, we must assume that the diagonalization subroutines produce eigenvectors \f$e^{i \theta_{\boldsymbol{k}}} \psi_{\boldsymbol{k}}\f$ with \f$\theta_{\boldsymbol{k}}\f$ a random number (even diagonalizing twice the Hamiltonian at the same wavevector may give different phases).

As a result, we patch the Quantum ESPRESSO code so that everytime the Hamiltonian is diagonalized in a non-self-consistent run, a gauge is fixed.
NOTA BENE: this implies that we need a nscf calculation before running ph.x!
We choose the gauge in this way.
Note that the wavefunction is expanded in a set of plane wave coefficients \f$ \psi_{\boldsymbol{k}} = \sum_{\boldsymbol{G}} c(\boldsymbol{G}) e^{i\boldsymbol{k}\cdot\boldsymbol{G}+i\boldsymbol{k}\cdot\boldsymbol{r}} \f$.
We fix the gauge such that the zero plane wave coefficient is real and positive: \f$ c(\boldsymbol{G}=0) = (C,0) \f$, with \f$C\f$ being a positive number.
While this gauge doesn't seem to have much physical value, is very simple to apply numerically.

There still is a catch in case of degenerate eigenvalues.
Let's make the example of a double-degenerate eigenvalue.
Let be \f$\psi_1\f$ and \f$\psi_2\f$ be two degenerate eigenvectors, then any \f$ \alpha \psi_1 + \beta \psi_2 \f$ is an eigenvector, as long as the coefficients are properly normalized.
More generally, in case of \f$n\f$ degenerate eigenvalues, with \f$\psi_i\f$ being the set of eigenvectors spanning the degenerate eigen-subspace, then we can rotate these wavefunctions as \f$ \theta_i = \sum_{j} R_{ij} \psi_{j} \f$, where \f$ R \f$ is any unitary matrix.
The unitary matrix is non-unique, and therefore it's much harder to fix a gauge for degenerate eigenvalues.

We therefore proceed by "fixing" the gauge in a way that is compatible with the wavefunction symmetries.
The wavefunction, at non-degenerate eigenvalues, satisfies two relations.
First, \f$ \psi_{\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}+\boldsymbol{G}}(\boldsymbol{r}) \f$ and then\f$ \psi_{S^{-1}\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(S\boldsymbol{r}) \f$, where \f$S\f$ is a symmetry operation of the crystal.
The symmetry operation \f$S\f$ consists of two operations \f$ S=\{R,t\} \f$: a rotation \f$R\f$ and a translation \f$\boldsymbol{t}\f$.
In terms of these quantities, the symmetry operation on the wavefunction is such that: \f$ \psi_{R\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(R^{-1}(\boldsymbol{r}-\boldsymbol{t})) \f$.
As a consequence, some further relations are:
\begin{equation}
\epsilon_{R\boldsymbol{k},n} = \epsilon_{\boldsymbol{k}n}
\end{equation}
\begin{equation}
\epsilon_{\boldsymbol{k}+\boldsymbol{K},n} = \epsilon_{\boldsymbol{k}n}
\end{equation}
\begin{equation}
c_{R\boldsymbol{k},n}(\boldsymbol{G}) = e^{-i(R\boldsymbol{k}+\boldsymbol{G}) \cdot \boldsymbol{t}} c_{\boldsymbol{k}n}(R^{-1}\boldsymbol{G})
\end{equation}
\begin{equation}
c_{\boldsymbol{k}+\boldsymbol{K},n}(\boldsymbol{G}) = c_{\boldsymbol{k}n}(\boldsymbol{G})
\end{equation}

Note, crucially, that the rotational/translational symmetries on the plane wave coefficients/ wavefunctions don't apply to degenerate states, due to the state-mixing problems mentioned above.
That is, the wavefunctions computed at rotated points might be mixed in different ways.

We therefore fix a gauge for degenerate states and restore the symmetries of the wavefunction by modifying the Quantum-ESPRESSO codebase with the following algorithm:

* The scf calculation is run using k-points in the irreducible wedge \f$ \{ k^{irr} \} \f$.
  After the diagonalization of the Hamiltonian,
  fix the gauge of non-degenerate eigenvectors setting c(G=0) to be real and positive.
  For degenerate eigenvalues, we only rotate c(G=0) for the first eigenvector in the degenerate subspace (whatever the value is).
  Save the wavefunction and its g-vectors (specifically, the arrays evc, igk_k).
  Save information on the symmetries that have been recognized.

* In any nscf/bands calculation, the Hamiltonian is diagonalized at a point k.
  If a non-degenerate band is found, rotate the wavefunctions with the c(G=0) plane wave coefficient as discussed above (both for non-degenerate bands and the first degenerate band).
  If a degenerate band is found:
  * Find the irreducible point, and the symmetry operation or translation \f$S\f$ such that \f$ Sk \in k^{irr} \f$.
  * Read from file the wavefunction at the irreducible k-point
  * Rotate the wavefunction. If S is the identity, this can be skipped.
  * Substitute the rotated wavefunction and energies in evc and et.
  * Recompute the ultrasoft/PAW part of the wavefunction (call calbec() )
    Otherwise, the augmented wavefunction may not use the correct phases.

Note that we are making the hypothesis that the nscf calculation is done using the same Monkhorst-Pack mesh that has been used in the scf run.

Note that both Wannier90 and ph.x must be run with the modified version of QE, so that the wavefunction gauge is fixed to be the same.

Note that the wavefunction, or g, doesn't need to be smooth with respect to \f$\boldsymbol{k}\f$: we would need to make sure that the derivative of the wavefunction with respect to \f$\boldsymbol{k}\f$ is continuous.

As a result, we fix both the gauge of the wavefunctions in degenerate subspaces and the gauge of the wavefunction across multiple \f$\boldsymbol{k}\f$ or \f$\boldsymbol{q}\f$ points.
To verify that the gauge has been fixed successfully, compare the electron-phonon coupling with and without symmetries, and they should be the same, see below.



\subsection COMPELPHQE2 Computational details, symmetries in Quantum ESPRESSO
Note: Abinit has a very well curated section on the symmetries of the wavefunction <https://docs.abinit.org/theory/wavefunctions/>

The phonon code can be used to compute the coupling \f$g(\boldsymbol{k},\boldsymbol{q})\f$, where k falls on a Monkhorst-Pack grid of points (nk1,nk2,nk3) and q falls on a Monkhorst-Pack mesh (nq1,nq2,nq3).
We require that both meshes are centered at the Gamma point, so that we have the Wannier90 matrices for the Bloch to Wannier rotation.
Given that the calculation is quite expensive, Quantum ESPRESSO uses symmetries to reduce the required amount of calculations.

As discussed above, we made sure that the set of wavefunctions obeys the relations: \f$ \psi_{R\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(R^{-1}(\boldsymbol{r}-\boldsymbol{t})) \f$.

Intuitively, the electron-phonon coupling itself should remain invariant under a symmetry operation: \f$ g(\boldsymbol{k},\boldsymbol{q}) = g(S\boldsymbol{k},S\boldsymbol{q}) \f$ and therefore, it should obey \f$g(S^{-1}\boldsymbol{k},\boldsymbol{q}) = g(\boldsymbol{k},S\boldsymbol{q})\f$. More systematically:

\f{eqnarray}{
g(\boldsymbol{k},S\boldsymbol{q})
&=& \big< \psi_{\boldsymbol{}k+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{S\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(S^{-1}\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(S\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(S\boldsymbol{r}) \big> \\
&=& \big< \psi_{S^{-1}\boldsymbol{k}+\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{S^{-1}\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& g(S^{-1}\boldsymbol{k},\boldsymbol{q})
\f}

Note two things: if the wavefunction doesn't rotate with the symmetries of the crystal (e.g. the gauge has not been fixed and degeneracies are not lifted), there will be phase factors hanging around, and the fourth equality in the expressions above doesn't hold.

Additionally, the translational invariance allows us to use the symmetry

\begin{equation} 
g(\boldsymbol{k},\boldsymbol{q}) = g(\boldsymbol{k}+\boldsymbol{G},\boldsymbol{q}+\boldsymbol{G}') \;,
\end{equation}

useful whenever a rotated point falls outside the Brillouin zone and must be folded back with an Umklapp vector \f$\boldsymbol{G}\f$.

The code ph.x uses two symmetries to reduce the list of \f$\boldsymbol{k}\f$ and \f$\boldsymbol{q}\f$ points.
First of all, ph.x only computes the coupling for the irreducible set of q wavevectors.
As a first guess, one may think that ph.x computes the coupling for all k points falling on a Monkhorst-Pack grid, for every irreducible q point.
However, at fixed irreducible \f$\boldsymbol{q}\f$ point, we don't need to compute all wavevectors \f$\boldsymbol{k}\f$.
In fact, consider a symmetry \f$S\f$ that sends the irreducible point \f$q\f$ to a reducible point \f$R\boldsymbol{q}=\boldsymbol{q}^*\f$ that are both on the Monkhorst-Pack mesh of q-points selected in input to ph.x.
While a wavevector \f$\boldsymbol{k}\f$ also falls on a Monkhors-Pack mesh, it may be that its rotation \f$R\boldsymbol{k}\f$ doesn't fall on the k-vector grid.
Therefore, we can discard the k-wavevectors of the grid that don't transform like \f$\boldsymbol{q}\f$ (for each irreducible q) and set their electron-phonon coupling to zero.
The ph.x code computes the coupling only for the pairs of \f$\boldsymbol{k}\f$ and \f$\boldsymbol{q}\f$ wavevectors that obey the same subset of symmetries, which can be rotated with the relations described above.
However, before testing this relation, we impose \f$\boldsymbol{k}\f$ to fall on a full grid.


\subsection PHSYMMS Computational details, phonon symmetries
We should not forget that also the phonon eigenvectors should satisfy the crystal symmetries when used for the Wannier transformation.
The symmetries of phonons are thoroughly discussed in this [reference] (https://link.aps.org/doi/10.1103/RevModPhys.40.1), from which we need just Eq. 2.33.
In detail, let the phonon eigenvector be \f$z_{\mu k j}(q)\f$, where \f$k\f$ is an atomic basis index, \f$\mu\f$ is a cartesian index, \f$q\f$ is the wavevector, and \f$j\f$ is the mode index.
If \f$S\f$ is a symmetry operation of the crystal, the phonon eigenvector rotates as:
\begin{equation}
\boldsymbol{q}' = S\boldsymbol{q}
\end{equation}
\begin{equation}
\omega_{j}(S\boldsymbol{q}) = \omega_{j}(\boldsymbol{q})
\end{equation}
\begin{equation}
z_{\mu K j}(S\boldsymbol{q}) = \sum_{\alpha} S_{\mu\nu} z_{\nu k j}(\boldsymbol{q}) \exp( i\boldsymbol{k} \cdot (S^{-1} R_{at}(K) - R_{at}(k)) )
\end{equation}
where \f$ R_{at} \f$ is the atomic position of an atom in the unit cell.
Furthermore, \f$K\f$ is the atomic basis index of the atom on which the atom \f$ k \f$ is transformed into upon the symmetry operation (since atoms of the same species are indistinguishable, they can be rotated into a different basis index, provided it's the same atomic species).
