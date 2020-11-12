@page thWannierElectrons Wannier interpolation methods

\section thElwHarmonic Wannier interpolation of band structure

A good review of the Wannier function formalism can be found [at this link] (https://arxiv.org/abs/1112.5411).

We assume that a third-party code has provided us with the single-particle hamiltonian in the wannier representation:

\begin{equation}
\langle \boldsymbol{0}n | H | \boldsymbol{R} m \rangle
\end{equation}
where \f$\boldsymbol{R}\f$ labels a bravais lattice vector in a supercell, \f$m\f$ and \f$n\f$ label two wannier functions.
A Wannier function is here denoted by the ket:
\begin{equation}
| \boldsymbol{R} m \rangle
\end{equation}

The Wannier interpolation requires first to transform to reciprocal space:
\begin{equation}
H_{\boldsymbol{k},nm}^W = \sum_{\boldsymbol{R}} e^{i \boldsymbol{k} \cdot \boldsymbol{R}} \langle \boldsymbol{0}n | H | \boldsymbol{R} m \rangle
\end{equation}
This sum is typically speeded-up summing over irreducible lattice vectors:
\begin{equation}
H_{\boldsymbol{k},nm}^W = \sum_{\boldsymbol{R}_{irr}} \frac{e^{i \boldsymbol{k} \cdot \boldsymbol{R}_{irr}} }{ d_{\boldsymbol{R}_{irr}}} \langle \boldsymbol{0} n | H | \boldsymbol{R} m \rangle
\end{equation}
where \f$d_{\boldsymbol{R}_{irr}}\f$ is the degree of degeneracy of the irreducible bravais lattice vector \f$\boldsymbol{R}_{irr}\f$.

This matrix is not diagonal, as we are working, typically, in the maximally localized wannier function gauge.
We thus diagonalize the matrix to pass to the Bloch representation:
\begin{equation}
H_{\boldsymbol{k},bb'}^B = [U_{\boldsymbol{k}}^\dagger H_{\boldsymbol{k}}^W U_{\boldsymbol{k}}]_{bb'} = \delta_{bb'} \epsilon_{\boldsymbol{k}b}
\end{equation}
Therefore finding the interpolated electronic energies at a generic \f$\boldsymbol{k}\f$ point.

The velocity operator is computed with the Hellman-Feynman theorem, e.g. along the x direction as:
\begin{equation}
v^x_{\boldsymbol{k}bb'} = [U_{\boldsymbol{k}}^\dagger \frac{H_{\boldsymbol{k}+\boldsymbol{\delta}_x}^W-H_{\boldsymbol{k}-\boldsymbol{\delta}_x}^W}{2 \delta_x} U_{\boldsymbol{k}}]_{bb'}
\end{equation}
Whenever we find a set of degenerate bands at a given k point, we diagonalize the velocity operator in the degenerate subset, in order to uniquely define the velocity operator.





@section thEpwInterpolation Wannier interpolation of electron-phonon coupling

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
\frac{1}{N_q}
\sum_{\boldsymbol{q}\nu} e^{-i\boldsymbol{q}\cdot\boldsymbol{R}_p} [u_{\boldsymbol{q}\kappa}^{\nu}]^{-1} \partial_{\boldsymbol{q}\nu} V(\boldsymbol{r})
\end{equation}


So, we first transform to Wannier space by:
\begin{equation}
g(\boldsymbol{R}_e,\boldsymbol{R}_p)
=
\frac{1}{N_q N_k}
\sum_{\boldsymbol{k}\boldsymbol{q}} e^{-i\boldsymbol{k}\cdot\boldsymbol{R}_e-i\boldsymbol{q}\cdot\boldsymbol{R}_p} U_{\boldsymbol{k}+\boldsymbol{q}}^\dagger g(\boldsymbol{k},\boldsymbol{q}) U_{\boldsymbol{k}} u_{\boldsymbol{q}}^{-1}
\end{equation}

Then, we interpolate to Bloch space

\begin{equation}
g(\boldsymbol{k},\boldsymbol{q})
=
\frac{1}{N_e}
\sum_{\boldsymbol{R}_e \boldsymbol{R}_p} e^{i\boldsymbol{k}\cdot\boldsymbol{R}_e+i\boldsymbol{q}\cdot\boldsymbol{R}_p} U_{\boldsymbol{k}+\boldsymbol{q}} g(\boldsymbol{R}_e,\boldsymbol{R}_p) U_{\boldsymbol{k}}^\dagger u_{\boldsymbol{q}}
\end{equation}


Details:
* the mesh of \f$\boldsymbol{k}\f$ and \f$\boldsymbol{q}\f$ points must be the same (or at least commensurate, so that we can map \f$\boldsymbol{k}+\boldsymbol{q}\f$ into the same \f$\boldsymbol{k}\f$ grid).






@subsubsection thCompW90 Computational details, Wannier90

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



\subsection thCOMPELPHQE Computational details, gauge fixing in Quantum ESPRESSO 

The interpolation procedure described above implicitely assumes that the wavefunction \f$\big|b\boldsymbol{k}\big>\f$ has a fixed gauge.
In fact, all the quantities above are complex numbers, and the wavefunction is defined within a phase (or, more generally, a unitary rotation).
In order for the interpolation to work, we must make sure that the wavefunction used for computing all the quantities above are exactly the same coefficient-wise, phase included, and make sure that pw.x, ph.x and wannier90.x operate on the same wavefunctions.

The problem comes from the arbitrariness of the phase choice of an eigenvector of a Hermitian matrix.
In details: let \f$H_{\boldsymbol{k}}\f$ be a Bloch Hamiltonian.
The DFT code will diagonalize the Hamiltonian and solvee \f$ H_{\boldsymbol{k}} \psi_{\boldsymbol{k}} = \epsilon_k \psi_{\boldsymbol{k}}\f$.
For each eigenvector \f$\psi_{\boldsymbol{k}}\f$, we can apply the transformation \f$\psi_{\boldsymbol{k}} \to e^{i \theta_{\boldsymbol{k}}} \psi_{\boldsymbol{k}}\f$ and still have \f$e^{i \theta_{\boldsymbol{k}}} \psi_{\boldsymbol{k}}\f$ an eigenvector.
Note also that the diagonalization may not have a strategy to fix the phase of the eigenvector: as a result, we may expect that every different run of a DFT code will generate a different phase, effectively behaving as a random number generator.

We thus patch the Quantum ESPRESSO code to fix a gauge of the wavefunction.
Additionally, we want to make sure that the wavefunction satisfies rotational symmetries, as this will help us reduce the number of calculations of the electron-phonon coupling at the DFT level.

In a plane-wave code, the wavefunction is expanded in a plane wave basis set as
\begin{equation}
\psi_{\boldsymbol{k}} = \sum_{\boldsymbol{G}} c(\boldsymbol{G}) e^{i\boldsymbol{k}\cdot\boldsymbol{G}+i\boldsymbol{k}\cdot\boldsymbol{r}}
\end{equation}
Quantum ESPRESSO, stores the plane wave coefficients in `evc(ig,ib)`, where `ib` is a band index and `ig` is a G-vector index.
Details are described in the source code, but keep in mind that `evc` is parallel-distributed over G-vectors, and that each k-point has a different order of G-vectors.
If we want to fix the gauge, we must operate on the plane wave coefficients.

The wavefunction satisfies some symmetries.
Let \f$S\f$ be a symmetry operation of the crystal.
A symmetry operation consists of a rotation \f$R\f$ and a fractional translation \f$t\f$, that leave the crystal invariant.
As the wavefunction must transform like the crystal, it can be shown that \f$ \psi_{R\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(R^{-1}(\boldsymbol{r}-\boldsymbol{t})) \f$.
From this symmetry property, one can verify that the following relations hold:
\begin{equation}
\epsilon_{R\boldsymbol{k},n} = \epsilon_{\boldsymbol{k}n}
\end{equation}
\begin{equation}
c_{R\boldsymbol{k},n}(\boldsymbol{G}) = e^{-i(R\boldsymbol{k}+\boldsymbol{G}) \cdot \boldsymbol{t}} c_{\boldsymbol{k}n}(R^{-1}\boldsymbol{G})
\end{equation}
Additionally, the wavefunction is periodic over the Brillouin zone, i.e. \f$ \psi_{k}(r) = \psi_{k+G'}(r) \f$.
From this, it follows that:
\begin{equation}
\epsilon_{\boldsymbol{k}+\boldsymbol{K},n} = \epsilon_{\boldsymbol{k}n}
\end{equation}
\begin{equation}
c_{\boldsymbol{k}+\boldsymbol{G}',n}(\boldsymbol{G})
=
c_{\boldsymbol{k}n}(\boldsymbol{G}+\boldsymbol{G}')
\end{equation}
Note: Abinit has a very well curated section on the <a href="https://docs.abinit.org/theory/wavefunctions/">symmetries of the wavefunction </a>

Before fixing the gauge, we also stress an additional problem: electronic degeneracy.
If two (or more) energy levels are degenerate, the wavefunction is only defined up to a unitary rotation.
In fact, let \f$i\f$ span the subspace of degenerate eigenvalues.
Then, the wavefunctions can be rotated as \f$ \tilde{\psi}_i = \sum_j U_{ij} \psi_j  \f$, with \f$ U \f$ any unitary matrix.
Therefore, when fixing the gauge, we must also deal with this problem: we must also mix the plane wave coefficients of different degenerate bands.

The algorithm to fix the gauge in Quantum ESPRESSO goes as follows:

* Run a scf calculation using the k-points in the irreducible wedge \f$ \{ k^{irr} \} \f$,
  setting the number of bands equal to what you want to use in both Wannier90 and ph.x.
  Right after the Hamiltonian is diagonalized at a given k-point (in file `PW/src/c_bands.f90`),
  and fix the gauge of non-degenerate eigenvectors by setting c(G=0) to be real and positive.
  For degenerate eigenvalues, set c(G=0)>0 only for the first band of the degenerate subspace.
  Save the wavefunction and its G-vectors (the arrays `evc`, `g_vectors`, and the mapping `igk_k`).
  
* During a ph.x calculation, or a nscf calculation before Wannier90, the codes ask to
  diagonalize the Hamiltonian at a point k (or k+q) that is commensurate with the grid of points
  used in the scf calculation.
  Given a point k, do:
  * find the irreducible point \f$ k^* \f$ that is symmetry-equivalent to the current point.
    If not found, block the code (the user has either messed symmetries or used wrong k/q meshes).
    Find also the symmetry operation S such that \f$ R k = k^* + K\f$,
    where \f$K\f$ is an Umklapp vector.
  * Read the wavefunction at \f$ k^* \f$.
  * Build `gmap`, a map between indices of two arrays of G-vectors such that
    \f$ G[i] = (R^{-1}G+K)[gmap(i)] \f$. This will help us apply the roto-translational symmetry.
  * Compute the roto-translated wavefunction \f$ \psi_{Rk} = \psi_{k^*+K} \f$
    using the relations on the plane-wave coefficients described above.
  
This would be enough, if the wavefunctions were exact.
Unfortunately, this procedure doesn't allow us to reconstruct the complete wavefunction.
In fact, the wavefunctions are typically expanded over a set of G-vectors such that \f$ |k+G|^2<E_{cut} \f$.
Therefore, the wavefunction can only be rotated for the intersecting set of G-vectors between the wavefunctions at the irreducible (reference) point and the roto-translated point.
We wouldn't have information for G-vectors outside this intersection and we would set them to zero, breaking the normalization condition.

We bypass this problem in this way.
Let \f$ \big| \psi^{QE} \big> \f$ be the wavefunction computed by QE at point k and \f$ \big| \psi^{rot} \big> \f$ the wavefunction we computed using the roto-translation of the irreducible point.
* Using the relation
  \f$ \big| \psi^{rot} \big>
  =
  \sum_{QE} \big< \psi^{rot} \big| \psi^{QE} \big>^* \big| \psi^{QE} \big>
  =
  U \big| \psi^{QE} \big>
  \f$
  to define a unitary matrix \f$ U \f$.
* On paper, \f$ U \f$ should be unitary, i.e. \f$ U U^{\dagger} = 1 \f$.
  But for the same problems of completeness of G-sphere, we have \f$ U U^{\dagger} = 1-\Delta \f$.
  With some manipulations,
  \f$ 1 = U U^{\dagger} + \Delta = U U^{\dagger} + U U^{\dagger} \Delta U U^{\dagger}
  = U ( 1 + U^{\dagger} \Delta U ) U^{\dagger}
  = U L L^{\dagger} U^{\dagger} \f$,
  where \f$L\f$ comes from the Cholesky decomposition of
  \f$( 1 + U^{\dagger} \Delta U ) = LL^{\dagger}\f$.
* Redefine \f$\tilde{U} = UL\f$ (this matrix is unitary by construction).
  Finally, the wavefunction at the point k is \f$ \tilde{U} \big| \psi^{QE} \big> \f$

This procedure has been implemented in QE, in the file `c_bands.f90`.

Note that there is a catch for entangled bands.
In building the unitary matrix \f$U\f$, we assumed completeness of the wavefunction set.
If you are Wannierizing disentangled bands, this is fine.
If you are trying to disentangle some bands, than it is possible that, by choosing the number of bands to be computed, we may cut through a group of degenerate bands.
If this happens, the last block of the matrix \f$U\f$ may not be unitary, not just because of numerical noise, but because of breaking the completeness relation.
We checked that, as long as you are discarding such bands in the disentangling procedure, the Wannierized wavefunctions should be fine.

Final comments:

1. In order to rotate the wavefunction, each MPI process needs to have enough memory to store
   the complete wavefunction (all G vectors) for a single band,
   i.e., each MPI process requires an additional \f$ 16 N_G \f$ Bytes of memory.

2. The lack of completeness implies that, as for any DFT calculation,
   one must converge the G-vectors cutoff (`ecutwfc` in QE).

3. The wavefunction, or g, even though it obeys symmetries,
   it isn't smooth with respect to \f$\boldsymbol{k}\f$.
   This is guaranteed by the maximally localized Wannier gauge
   (which in the reciprocal space guarantees continuity with respect to k).

4. Currently we don't support spin, but we will add it soon (must include a few more symmetries).






\subsection thCOMPELPHQE2 Computational details, symmetries in Quantum ESPRESSO

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



\section thElScatt Electron BTE

Let \f$ f_{\nu} \f$ be the out-of-equilibrium electron occupation number, where \f$ \nu = (\boldsymbol{k},b) \f$ labels both electronic wavevectors and band index (i.e. the single-particle Bloch numbers).
First, we rewrite the occupation number as:

\f[
f_{\lambda} = \bar{f}_{\lambda} + \bar{f}_{\lambda}(1-\bar{f}_{\lambda}) \delta f_{\lambda} 
\f]
where \f$ \bar{f}_{\lambda} \f$ is the Fermi--Dirac distribution function and we introduced \f$ \delta f_{\lambda} \f$ as the canonical distribution function.

The linearized electronic BTE can be written as

\f[
 - e \boldsymbol{v}_{\lambda} \cdot \boldsymbol{E} \frac{\partial \bar{f}_{\lambda}}{\partial \epsilon} + \boldsymbol{v}_{\lambda} \cdot \boldsymbol{\nabla} T \frac{\partial \bar{f}_{\lambda}}{\partial T} =
- \sum_{\lambda'} A_{\lambda\lambda'} \delta f_{\lambda'}
\f]
where the first term describes the diffusion due to an externally applied electric field \f$ \boldsymbol{E} \f$, the second  the diffusion due to a temperature gradient, and the third term is the linearized scattering operator.

The scattering matrix \f$ A_{\lambda,\lambda'} \f$ can be computed as
\f[
A_{\boldsymbol{k}b,\boldsymbol{k}'b'} = \frac{1}{V N_k} \sum_{s, \boldsymbol{q}}
2 \pi
|g_{bb'\nu}(\boldsymbol{k},\boldsymbol{k}')|^2
\times
\bigg[
\bar{f}_{\boldsymbol{k}b}(1-\bar{f}_{\boldsymbol{k}'b'}) \bar{n}_{\boldsymbol{q}\nu}
\delta(\epsilon_{\boldsymbol{k}b} + \hbar \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}'b'})
+
\bar{f}_{\boldsymbol{k}'b'}(1-\bar{f}_{\boldsymbol{k}b}) \bar{n}_{\boldsymbol{q}\nu}
\delta(\epsilon_{\boldsymbol{k}b} - \hbar \omega_{\boldsymbol{q}\nu} - \epsilon_{\boldsymbol{k}'b'})
\bigg]
\delta(\boldsymbol{k}-\boldsymbol{k}'+\boldsymbol{q})
\f]

This quantity can be computed knowing all the interpolation techniques on phonon energies, electronic energies and the electron-phonon coupling.
Please note that, for convenience, here we use a coupling defined as
\f[
g_{bb'\nu}(\boldsymbol{k},\boldsymbol{k}')
=
g_{b'b\nu}(\boldsymbol{k},\boldsymbol{q})
\f]
where the latter can be interpolated as described above.
The Dirac-delta conserving momentum is enforced exactly, since we are using points on a uniform grid centered at gamma.
The Dirac-delta conserving energy is instead with a Gaussian function, as described in the section @ref thSMEARING for the phonon BTE.





\section thOnsager Onsager coefficients
We make the hypothesis that the response is linear in the external fields:
\f[
\delta f_{\lambda} = \sum_{i} \delta^i f^E_{\lambda} E_i + \delta^i f^T_{\lambda} \nabla_i T 
\f]

After computing the out-of-equilibrium population, the charge and heat flux density can be computed as:
\f[
\boldsymbol{J} = \frac{e g_s}{V N_k} \sum_{\lambda} \boldsymbol{v}_{\lambda} f_{\lambda}
\f]
and
\f[
\boldsymbol{Q} = \frac{g_s}{V N_k} \sum_{\lambda} (\epsilon_{\lambda}-\mu) \boldsymbol{v}_{\lambda} f_{\lambda}
\f]
where \f$ g_s\f$ is the spin degeneracy.

Thanks to the decomposition, we can write
\f[
\boldsymbol{J} = L_{EE} \boldsymbol{E} + L_{ET} \boldsymbol{\nabla} T
\f]
\f[
\boldsymbol{Q} = L_{TE} \boldsymbol{E} + L_{TT} \boldsymbol{\nabla} T
\f]

The electrical conductivity \f$ \sigma \f$, the thermal conductivity \f$ k \f$, the Seebeck coefficient \f$ S \f$ and the mobility \f$ \mu \f$ are:
\f[
\sigma = L_{EE}
\f]
\f[
k = L_{TT} - L_{TE} L_{EE}^{-1} L_{ET}
\f]
\f[
S = - L_{EE}^{-1} L_{ET}
\f]
\f[
\mu = \frac{\sigma}{d}
\f]
where \f$ d\f$ is the carriers'/doping concentration.




\section thElRTA Electron transport in relaxation time approximation
At this simple level of theory, we define the electron lifetime as:
\f[
A_{ \boldsymbol{k}b,\boldsymbol{k}b } = \frac{\bar{f}_{\boldsymbol{k}b}(1-\bar{f}_{\boldsymbol{k}b})}{ \tau_{\boldsymbol{k}b} }
\f]
Next, we approximate the scattering matrix as diagonal, so that the BTE becomes:

\f[
 - e \boldsymbol{v}_{\lambda} \cdot \boldsymbol{E} \frac{\partial \bar{f}_{\lambda}}{\partial \epsilon} + \boldsymbol{v}_{\lambda} \cdot \boldsymbol{\nabla} T \frac{\partial \bar{f}_{\lambda}}{\partial T} =
- \frac{\bar{f}_{\lambda}(1-\bar{f}_{\lambda})}{ \tau_{\lambda} } \delta f_{\lambda}
\f]

Solving separately the response to the electric field and the thermal gradient, we find
\f[
\delta^i f^E_{\lambda} = - e v^i_{\lambda} \frac{1}{k_B T} \tau_{\lambda}
\f]
\f[
\delta^i f^T_{\lambda} = v^i_{\lambda} \frac{\epsilon_{\lambda}}{k_B T^2} \tau_{\lambda}
\f]



@section thElIterative Iterative solution to the electron BTE - Omini Sparavigna method

This is an adaptation of the Omini-Sparavigna method to electrons.
<blockquote>
Generally, we recommend the variational method over this. 
</blockquote>
To better understand this method, please have a look first at the phonon counterpart @ref thPHITER.

The BTE consists in two linear algebra problems:
\f[
m^{i}_{\lambda} = - \sum_{\lambda'} A_{\lambda\lambda'} \delta f_{\lambda}^E
\f]
\f[
n^{i}_{\lambda} = - \sum_{\lambda'} A_{\lambda\lambda'} \delta f_{\lambda}^T
\f]
where
\f[
m^{i}_{\lambda} = - e v_{\lambda}^i \frac{\partial \bar{f}_{\lambda}}{\partial \epsilon}
\f]
\f[
n^{i}_{\lambda} = v_{\lambda}^i \frac{\partial \bar{f}_{\lambda}}{\partial T}
\f]

The iterative scheme consists in solving iteratively this two independent linear algebra problems with geometric series:
\f[
\delta^i f^E_{K} = \sum_{K} \left(-\frac{1}{\boldsymbol{A}^{\mathrm{out}}}  \boldsymbol{A}^{\mathrm{in}}\right)^{K} \frac{1}{\boldsymbol{A}^{\mathrm{out}}} \:  m^i
\f]
and
\f[
\delta^i f^T_K = \sum_{K} \left(-\frac{1}{\boldsymbol{A}^{\mathrm{out}}}  \boldsymbol{A}^{\mathrm{in}}\right)^{K} \frac{1}{\boldsymbol{A}^{\mathrm{out}}} \:  n^i
\f]
where \f$ K \f$ is an iteration index, \f$ A^{in} \f$ is the off-diagonal part of the scattering matrix, and \f$ A^{out} \f$ is the diagonal part of the scattering matrix.
Note that, like any geometric series, this algorithm may not converge.
In the code, the two problems are solved together, as we compute the action on the two different vectors at the same time.






@section thElVariational Variational solution to the electron BTE

As seen above, the electron solvers to the BTE are identical to the phonon case.
The only difference is that we need to solve two problems simultaneously, one for the electric field response and one for the response to the thermal gradient.

For the variational method, we can define the variational thermal conductivity, in closed-circuit conditions, as:
\f[
k^\mathrm{V}(\delta f^T) = - 2 \mathcal{T}({\delta f^T})
\f]
where
\f[
\mathcal{T}(\delta f^T) = \frac{1}{2} \sum_{\lambda \lambda'} {\delta f^T_{\lambda}} \cdot{\boldsymbol A_{\lambda\lambda'}} {\delta f^T_{\lambda'}} - \sum_{\lambda} {\boldsymbol n_{\lambda}} \cdot {\delta f^T_{\lambda}}
\f]

The variational electrical conductivity is defined similarly as:
\f[
\sigma^\mathrm{V}(\delta f^E) = 2 \mathcal{E}({\delta f^E})
\f]
where
\f[
\mathcal{E}(\delta f^E) = \frac{1}{2} \sum_{\lambda \lambda'} {\delta f^E_{\lambda}} \cdot{\boldsymbol A_{\lambda\lambda'}} {\delta f^E_{\lambda'}} - \sum_{\lambda} {\boldsymbol m_{\lambda}} \cdot {\delta f^E_{\lambda}}
\f]


These two functionals are the minimization targets of a conjugate gradient method.
Knowing this, the variational method is exactly the same as the phonon case described in section @ref thPHITER, with the proper substitution of the vector \f$b\f$ with either \f$m\f$ or \f$n\f$.

As in the case of the Omini-Sparavigna method, we solve the two equations (response to electric field and thermal gradient) at the same time, as it allows us to minimize the number of times the scattering matrix is evaluated (the most expensive step).








@section thElRelaxons Relaxons solution to the electronic BTE
In this scheme, we use an algebraic solution to the BTE, solving the equation in the eigenvector basis.
We first diagonalize the scattering matrix:
\begin{equation}
\frac{1}{N_k} \sum_{\lambda'} A_{\lambda\lambda'} \theta_{\lambda'\alpha} = \frac{1}{\tau_{\alpha}} \theta_{\lambda\alpha}
\end{equation}
where \f$ \theta \f$ are eigenvectors, \f$ \alpha \f$ are eigenvalue indices, and \f$ \frac{1}{\tau_{\alpha}} \f$ are eigenvalues.
We first build the auxiliary quantities:
\f[
\delta^i f^E_{\alpha} = - \sum_{\lambda} \frac{\partial \bar{f}_{\lambda}}{\partial \epsilon} v_{\lambda}^i  \theta_{\lambda \alpha} \tau_{\alpha}
\f]
\f[
\delta^i f^T_{\alpha} = \sum_{\lambda} \frac{\partial \bar{f}_{\lambda}}{\partial T} v_{\lambda}^i  \theta_{\lambda \alpha} \tau_{\alpha}
\f]
From these, we can compute the solutions of the BTE as:
\f[
\delta f^E_{\lambda} = \frac{1}{N_k V} \sum_{\alpha} f^E_{\alpha} \theta_{\lambda \alpha}
\f]
\f[
\delta f^T_{\lambda} = \frac{1}{N_k V} \sum_{\alpha} f^T_{\alpha} \theta_{\lambda \alpha}
\f]





\section thElWigner Wigner correction to the BTE
The Wigner transport equation is
\f[
\frac{\partial f_{bb'}(\boldsymbol{x},\boldsymbol{k},t)}{\partial t}
+
 \frac{i}{\hbar} \Big[ \mathcal{E}(\boldsymbol{k}) + \boldsymbol{D}(\boldsymbol{k})\cdot\boldsymbol{E} , f(\boldsymbol{x},\boldsymbol{k},t) \Big]_{bb'}
+
 \frac{1}{2} \Big\{ \boldsymbol{v}(\boldsymbol{k}) , \cdot \frac{\partial f(\boldsymbol{x},\boldsymbol{k},t)}{\partial \boldsymbol{x}} \Big \}_{bb'}
-
 e \boldsymbol{E} \cdot \frac{\partial f_{bb'}(\boldsymbol{x},\boldsymbol{k},t)}{\partial \boldsymbol{k}}
=
-\frac{\partial f_{bb'}(\boldsymbol{x},\boldsymbol{k},t)}{\partial t} \bigg|_{coll} 
\f]
where \f$ f_{bb'}(\boldsymbol{x},\boldsymbol{k},t) \f$ is the Wigner distribution function, \f$ \{ \cdot,\cdot \} \f$ indicates an anticommutator, \f$ [ \cdot,\cdot ] \f$ indicates a commutator, \f$ v_{bb'}(\boldsymbol{k}) \f$ is the velocity operator, and we defined the matrix \f$ \mathcal{E}(\boldsymbol{k})_{bb'} = \delta_{bb'} \epsilon_{\boldsymbol{k}b} \f$ and \f$ \mathcal{D}(\boldsymbol{k})_{bb'} = (1-\delta_{bb'}) d_{\boldsymbol{k}bb'} \f$ is a matrix of electronic dipoles.
The electronic dipole can be computed as:
\f[
\boldsymbol{d}_{\boldsymbol{k},bb'}
=
- i e \frac{\boldsymbol{v}_{bb'}(\boldsymbol{k})}{\epsilon_{b}(\boldsymbol{k})-\epsilon_{b'}(\boldsymbol{k})}  , \quad \text{for }b \neq b'
\f]

The scattering operator acts on the diagonal Wigner distribution as the BTE scattering operator, instead it acts on the off-diagonal components with a decay term:
\f[
\frac{\partial f_{bb'}(\boldsymbol{x},\boldsymbol{k},t)}{\partial t} \bigg|_{coll} 
=
(1-\delta_{bb'}) \frac{\Gamma_{b}(\boldsymbol{k}) + \Gamma_{b'}(\boldsymbol{k})}{2} f_{bb'}(\boldsymbol{x},\boldsymbol{k},t)
+
\delta_{bb'} \frac{1}{V}
\sum_{\boldsymbol{k}'b'} A_{\boldsymbol{k}b,\boldsymbol{k}'b'} f_{b'b'}(\boldsymbol{x},\boldsymbol{k}',t)
\f]
where \f$ \Gamma_b(\boldsymbol{k}) = \frac{2\pi}{\tau_{\boldsymbol{k}b}} \f$ are the electronic linewidths.

To solve the Wigner transport equation, just like we did for the BTE, we assume linear response and separate the response to electric field and thermal gradient \f$ f = f^E E + f^T \nabla T \f$.
The diagonal part of the Wigner transport equation is exactly equal to the BTE, and can be solved using one of solvers described above.
The off-diagonal part of the Wigner distribution function can be solved easily with a little algebraic manipulation.

The transport coefficients are defined as:
\f[
L_{EE}^{ij} = 
\frac{e g_s}{V N_k} \sum_{\boldsymbol{k}b} \frac{1}{2} \Big\{ v^i(\boldsymbol{k}) , f^{E_j}(\boldsymbol{k}) \Big\}_{bb}
\f]
\f[
L_{ET}^{ij} = 
\frac{e g_s}{V N_k} \sum_{\boldsymbol{k}b} \frac{1}{2} \Big\{ v^i(\boldsymbol{k}) , f^{T_j}(\boldsymbol{k}) \Big\}_{bb}
\f]
\f[
L_{TE}^{ij} = 
\frac{g_s}{V N_k} 
\sum_{\boldsymbol{k}b} 
\big( \epsilon_{b}(\boldsymbol{k})-\mu \big)
\frac{1}{2} \Big\{ v^i(\boldsymbol{k}) , f^{E_j}(\boldsymbol{k}) \Big\}_{bb}
\f]
\f[
L_{TT}^{ij} = 
\frac{g_s}{V N_k} 
\sum_{\boldsymbol{k}b} 
\big( \epsilon_{b}(\boldsymbol{k})-\mu \big)
\frac{1}{2} \Big\{ v^i(\boldsymbol{k}) , f^{T_j}(\boldsymbol{k}) \Big\} _{bb}
\f]
