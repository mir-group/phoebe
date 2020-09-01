@page Theory Theory
@section ELPHC Electron-phonon coupling, Wannier interpolation

This is a summary of the strategy for interpolation adopted by Phoebe, as detailed also at this [reference] (https://journals.aps.org/prb/abstract/10.1103/PhysRevB.76.165108).

The coupling to be used for the calculation of scattering rates is:


\begin{equation}
g^{SE}_{b'b,\nu} (k,q) = \bigg( \frac{1}{2 m \omega_{q\nu}} \bigg)^{1/2} g_{b'b,\nu} (k,q)
\end{equation}

What we need to interpolate instead is:
\begin{equation}
g_{b'b,\nu} (k,q) = \big<b'k+q \big| \partial_{q\nu}V \big| bk \big>
\end{equation}

First recall the relation between the wavefunction in Wannier and Bloch representation.

\begin{equation}
\big|mR_e\big> = \sum_{bk} e^{-ikR_e} U_{mb,k} \big|bk\big>
\end{equation}

\begin{equation}
\big|bk\big> = \frac{1}{N_e} \sum_{mR_e} e^{ikR_e} U_{bm,k}^\dagger \big|mR_e\big>
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


@subsection COMPELPH Computational details, Wannier90

First, the mesh of k and q points must be commensurate.
In fact, the matrices \f$U\f$ are computed on a grid of \f$k\f$ points.
It is thus necessary that \f$k+q\f$ falls on the same grid of points.


To be precise, note that we use two different sets of \f$U\f$ matrices.
When interpolating from the Wannier to the Bloch space, we use the U matrices computed at an arbitrary grid of points obtained by diagonalizing the Hamiltonian matrix in the Wannier space.
When first transforming from the Bloch space (the coupling with an ab-initio code) to the Wannier space, we use the \f$U\f$ matrices computed by Wannier90.
This is necessary, because the coupling \f$g\f$ is a complex function, and we must rotate it using the same wavefunction gauge used to build the maximally localized wannier functions.
Moreover, Wannier90 may use the disentanglement procedure.
In that case, the Bloch to Wannier transformation is:
\begin{equation}
\big|mR_e\big> = \sum_{k b b'} e^{-ikR_e} U^{rot}_{mb',k} U^{dis}_{b'b,k} \big|bk\big>
\end{equation}
where the number of disentangled bands \f$b'\f$ is smaller than the number of entangled bands \f$b\f$.
Therefore, we rotate the electron-phonon coupling from the Bloch to Wannier space using the entangled number of bands.
Wannier90 prints the two different \f$U\f$ matrices, and one can just multiply them to get the transformation matrix.

As a further minor detail, remember that some bands (like deep core bands) may be excluded from the Wannierization procedure (through the keyword exclude-indeces), so that there may be an offset in the band index of U and g.



\subsection COMPELPHQE Computational details, gauge in Quantum ESPRESSO 

The interpolation procedure described above implicitely assumes that the wavefunction \f$\big|bk\big>\f$ has a fixed gauge.
This is because the rotation \f$U\f$ must act on the same wavefunction used to compute the coupling \f$g\f$.

If no gauge is imposed on the Bloch wavefunction, \f$g\f$ becomes a white noise.
In details: consider the Bloch Hamiltonian \f$H_k\f$.
A DFT code using periodic boundary conditions will diagonalize it to find \f$ H_k \psi_k = \epsilon_k \psi_k\f$.
However, we can adopt the transformation \f$\psi_k \to e^{i \theta_k} \psi_k\f$ and still have \f$e^{i \theta_k} \psi_k\f$ be an eigenvector.
Specifically, the numerical diagonalization procedure doesn't by itself fix a gauge: in fact, we must assume that the diagonalization subroutines produce eigenvectors \f$e^{i \theta_k} \psi_k\f$ with \f$\theta_k\f$ a random number (even diagonalizing twice the Hamiltonian at the same wavevector may give different phases).

As a result, we patch the Quantum ESPRESSO code so that everytime the Hamiltonian is diagonalized in a non-self-consistent run, a gauge is fixed.
NOTA BENE: this implies that we need a nscf calculation before running ph.x!
We choose the gauge in this way.
Note that the wavefunction is expanded in a set of plane wave coefficients \f$ \psi_k = \sum_G c(G) e^{ikG+ikr} \f$.
We fix the gauge such that the largest plane wave coefficient is real: \f$ max( |c(G)| ) = (C,0) \f$, with \f$C\f$ being a positive number.
While this gauge doesn't seem to have much physical value, is very simple to apply.

TODO: there is an ambiguity in how this max is found, in case there are degenerate maxs. This ambiguity might come into play if we change the parallelization between different runs. A similar problem happens below which may appear with different parallelization schemes in different runs.

There still is a catch in case of degenerate eigenvalues.
Let's make the example of a double-degenerate eigenvalue.
Let be \f$\psi_1\f$ and \f$\psi_2\f$ be two degenerate eigenvectors, then any \f$ \alpha \psi_1 + \beta \psi_2 \f$ is an eigenvector.
For the same arguments as above, the diagonalization procedure essentially assigns a random phase to the every degenerate subspace.
Therefore, we need to fix this degree of freedom.
To this aim, after the gauge on c(G) is imposed, we look for degenerate subspaces.
In each degenerate subspace of size \f$D\f$, we compute the DxD matrix \f$ F = \sum_{i=1,D; j=1,D} c(G,i) K(i,j) c(G',j) \f$ where we only use the first \f$D\f$ plane wave coefficients found on the root node, i and j are band indeces in the degenerate subspace, and \f$ K \f$ is a Hermitian DxD matrix with pseudo-random values (pseudo-random so that we can systematically recreate it).
We diagonalize the \f$F\f$ matrix and check that the eigenvalues of this new matrix are non-degenerate.
Finally, we rotate the original degenerate wavefunctions \f$\psi_i\f$ using the matrix of eigenvectors of \f$F\f$.
This procedure should be equivalent to lifting the degeneracy of the Hamiltonian by applying a small perturbation, but more computationally efficient.

As a result, we fix both the gauge of the wavefunctions in degenerate subspaces and the gauge of the wavefunction across multiple k or q points.
Now, the \f$g\f$ coupling, computed at fixed gauge, is a smooth function of k and q, which can be interpolated.
Note that also Wannier90 must be used starting from the gauge-fixed wavefunctions.



\subsection COMPELPHQE2 Computational details, symmetries in Quantum ESPRESSO
The phonon code can be used to compute the coupling \f$g(k,q)\f$, where k falls on a Monkhorst-Pack grid of points (nk1,nk2,nk3) and q falls on a Monkhorst-Pack mesh (nq1,nq2,nq3).
We require that both meshes are centered at the Gamma point, so that we have the Wannier90 matrices for the Bloch to Wannier rotation.
Given that the calculation is quite expensive, Quantum ESPRESSO uses symmetries to reduce the required amount of calculations.

After the wavefunction gauge is fixed (and only after) we have two nice symmetries: \f$ \psi_k(r) = \psi_{k+G}(r) \f$ and \f$ \psi_{S^{-1}k}(r) = \psi_{k+G}(Sr) \f$, where \f$S\f$ is a symmetry operation of the crystal.
The electron-phonon coupling itself should remain invariant under a symmetry operation: \f$ g(k,q) = g(Sk,Sq) \f$. Therefore, it follows that

\begin{equation}
g(S^{-1}k,q) = g(k,Sq) \;.
\end{equation}

Additionally, the Bloch condition allows us to use the symmetry

\begin{equation} 
g(k,q) = g(k+G,q+G') \;,
\end{equation}

useful whenever a rotated point falls outside the Brillouin zone and must be folded back with an Umklapp vector \f$G\f$.

The code ph.x uses two symmetries to reduce the list of k and q points.
First of all, ph.x only computes the coupling for the irreducible set of q wavevectors.
As a first guess, one may think that ph.x computes the coupling for all k points falling on a Monkhorst-Pack grid, for every irreducible q point.
However, at fixed irreducible q point, we don't need to compute all wavevectors k.
In fact, consider a symmetry \f$S\f$ that sends the irreducible point \f$q\f$ to a reducible point \f$Sq=q^*\f$ that are both on the Monkhorst-Pack mesh of q-points selected in input to ph.x.
While a wavevector k also falls on a Monkhors-Pack mesh, it may be that its rotation \f$Sk\f$ doesn't fall on the k-vector grid.
Therefore, we can discard the k-wavevectors of the grid that don't transform like q (for each irreducible q) and set their electron-phonon coupling to zero.
The ph.x code computes the coupling only for the pairs of k and q wavevectors that obey the same subset of symmetries, which can be rotated with the relations described above.
However, before testing this relation, we impose k to fall on a full grid.