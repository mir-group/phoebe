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



\subsection COMPELPHQE Computational details, gauge in Quantum ESPRESSO 

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
We fix the gauge such that the first plane wave coefficient is real: \f$ c(\boldsymbol{G}=0) = (C,0) \f$, with \f$C\f$ being a positive number.
While this gauge doesn't seem to have much physical value, is very simple to apply numerically.

There still is a catch in case of degenerate eigenvalues.
Let's make the example of a double-degenerate eigenvalue.
Let be \f$\psi_1\f$ and \f$\psi_2\f$ be two degenerate eigenvectors, then any \f$ \alpha \psi_1 + \beta \psi_2 \f$ is an eigenvector.
For the same arguments as above, the diagonalization procedure essentially assigns a random phase to the every degenerate subspace.
Therefore, we need to fix this degree of freedom.
To this aim, we must define a perturbation matrix, compute the expectation values against the wavefunction, diagonalize it, and use its eigenvectors to rotate the degenerate wavefunctions.
To be more computationally efficient, the wavefunction is restricted to the \f$\boldsymbol{G}=0\f$ plane wave coefficient.
Working on \f$\boldsymbol{G}=0\f$ also allows us to ignore the rotations on G-vectors.
The perturbation matrix must be built in a pseudo-random fashion (so that it's always the same for a fixed degenerate dimension).
The perturbation matrix must also reflect the crystal symmetry, so that the degeneracy of symmetry-equivalent k points is lifted in the same way.
For this, we take inspiration from the Fourier interpolation and multiply the constant perturbation matrix with a function \f$ S(\boldsymbol{k}) = \frac{1}{n_{\boldsymbol{R}}} \sum_{\boldsymbol{R}} e^{i \boldsymbol{k} \cdot \boldsymbol{R}}  \f$ is a function with the same symmetries of the crystal, with \f$\boldsymbol{R}\f$ being a set of irreducible Bravais lattice vectors (of an arbitrarily chosen supercell) and \f$n_R\f$ is the degeneracy of the Bravais lattice vectors.

Summarising, the algorithm goes as follow (we first lift degeneracies and then fix gauge):
* Find, in Quantum ESPRESSO, the place right after the Hamiltonian \f$H_{\boldsymbol{k}}\f$ has been diagonalized at each wavevector, and we have both the set of eigenvalues and eigenvectors \f$ c_{\boldsymbol{k},b}(\boldsymbol{G}) \f$ and \f$ \epsilon_{\boldsymbol{k},b} \f$.
* Look for degenerate subspaces
* For each degenerate subspace (of space size \f$d\f$), pick the set of coefficients \f$ c = c_{\boldsymbol{k},i}(\boldsymbol{G}=0) \f$, for i in the degenerate subspace.
* Build a perturbation matrix \f$ \Delta H_{\boldsymbol{k}} = F S(\boldsymbol{k}) \f$, where \f$ F \f$ is a pseudo-random \f$ d\times d\f$ hermitian matrix with non-degenerate eigenvalues, and \f$ S(\boldsymbol{k}) = \frac{1}{n_R} \sum_{\boldsymbol{R}} e^{i \boldsymbol{k} \cdot \boldsymbol{R}}  \f$ is a function with the same symmetries of the crystal, with \f$\boldsymbol{R}\f$ being a set of irreducible Bravais lattice vectors (of an arbitrarily chosen supercell) and \f$n_R\f$ is the degeneracy of the Bravais lattice vectors.
* Build the \f$ d \times d \f$ matrix \f$ D = ( c^* \otimes c ) \cdot \Delta H \f$, where \f$ \otimes \f$ is the tensor product.
* Diagonalize the matrix D, \f$ D\theta=\lambda\theta \f$, check that the eigenvalues have lost their degeneracy, and use the eigenvectors to rotate the plane wave coefficients at \f$ \boldsymbol{G}=0 \f$, i.e. \f$ c^{rot} = \theta^{\dagger} c \f$.
* After all degeneracies have been lifted, consider again the list of \f$\boldsymbol{G}=0\f$ coefficients \f$ c = c_{\boldsymbol{k},b}(\boldsymbol{G}=0)\f$.
For each band, find the phase that makes the list of \f$c\f$ coefficients real and positive.
* (Broadcast these phases, and) rotate the full list of plane wave coefficients.

Note that only root node needs to work, and only broadcasts the rotations twice (once for the degeneracy lifting, and once for the rotation. Considering that the degenerate spaces are typically small (i.e. typically dimensions smaller than 10), and that we can limit the Bravais lattice vectors to be less than 100, the whole procedure is negligibly expensive compared to the diagonalization of the \f$ H_{\boldsymbol{k}}(\boldsymbol{G},\boldsymbol{G}')\f$ Hamiltonian.

Note that Wannier90 must be run starting from the gauge-fixed wavefunctions.

Note that the wavefunction, or g, doesn't need to be smooth with respect to \f$\boldsymbol{k}\f$: we would need to make sure that the derivative of the wavefunction with respect to \f$\boldsymbol{k}\f$ is continuous.

As a result, we fix both the gauge of the wavefunctions in degenerate subspaces and the gauge of the wavefunction across multiple \f$\boldsymbol{k}\f$ or \f$\boldsymbol{q}\f$ points.
To verify that the gauge has been fixed successfully, compare the electron-phonon coupling with and without symmetries, and they should be the same, see below.



\subsection COMPELPHQE2 Computational details, symmetries in Quantum ESPRESSO
Note: Abinit has a very well curated section on the symmetries of the wavefunction <https://docs.abinit.org/theory/wavefunctions/>

The phonon code can be used to compute the coupling \f$g(\boldsymbol{k},\boldsymbol{q})\f$, where k falls on a Monkhorst-Pack grid of points (nk1,nk2,nk3) and q falls on a Monkhorst-Pack mesh (nq1,nq2,nq3).
We require that both meshes are centered at the Gamma point, so that we have the Wannier90 matrices for the Bloch to Wannier rotation.
Given that the calculation is quite expensive, Quantum ESPRESSO uses symmetries to reduce the required amount of calculations.

After the wavefunction gauge is fixed (and only after) we have two nice symmetries: \f$ \psi_{\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}+\boldsymbol{G}}(\boldsymbol{r}) \f$ and \f$ \psi_{S^{-1}\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(S\boldsymbol{r}) \f$, where \f$S\f$ is a symmetry operation of the crystal.
The symmetry operation \f$S\f$ consists of two operations: a rotation \f$R\f$ and a translation \f$\boldsymbol{t}\f$.
In terms of these quantities, the symmetry operation on the wavefunction is such that: \f$ \psi_{R\boldsymbol{k}}(\boldsymbol{r}) = \psi_{\boldsymbol{k}}(R^{-1}(\boldsymbol{r}-\boldsymbol{t})) \f$.
The most important relations are:
\begin{equation}
\epsilon_{R\boldsymbol{k}} = \epsilon_{\boldsymbol{k}}
\end{equation}
\begin{equation}
c_{R\boldsymbol{k}}(\boldsymbol{G}) = e^{-i(R\boldsymbol{k}+\boldsymbol{G}) \cdot \boldsymbol{t}} c_{\boldsymbol{k}}(R^{-1}\boldsymbol{G})
\end{equation}


Assume that the wavefunction gauge obeys the symmetries of the crystal and that degeneracies have been removed.
Intuitively, the electron-phonon coupling itself should remain invariant under a symmetry operation: \f$ g(\boldsymbol{k},\boldsymbol{q}) = g(S\boldsymbol{k},S\boldsymbol{q}) \f$ and therefore, it should obey \f$g(S^{-1}\boldsymbol{k},\boldsymbol{q}) = g(\boldsymbol{k},S\boldsymbol{q})\f$. More systematically:

\f{eqnarray}{
g(\boldsymbol{k},S\boldsymbol{q})
&=& \big< \psi_{\boldsymbol{}k+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{S\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(S^{-1}\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(S\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(S\boldsymbol{r}) \big> \\
&=& \big< \psi_{S^{-1}\boldsymbol{k}+\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{S^{-1}\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& g(S^{-1}\boldsymbol{k},\boldsymbol{q})   \quad \quad \text{(Wrong! See below)}
\f}

Note two things: if the wavefunction doesn't rotate with the symmetries of the crystal (e.g. the gauge has not been fixed and degeneracies are not lifted), there will be phase factors hanging around, and the fourth equality in the expressions above doesn't hold.

Now, note that the wavefunction, after fixing the gauge as described in the previous paragraph, doesn't obey the symmetries.
In fact, since we imposed that the plane wave coefficient at \f$\boldsymbol{G}=0\f$ is real and  positive.
Therefore, while we would like to have \f$c_{R\boldsymbol{k}}(\boldsymbol{G}=0) = e^{-iR\boldsymbol{k}\cdot\boldsymbol{t}} c_{\boldsymbol{k}}(\boldsymbol{G}=0)\f$ , we have instead imposed the rule \f$c_{R\boldsymbol{k}}(\boldsymbol{G}=0) = c_{\boldsymbol{k}}(\boldsymbol{G}=0)\f$.
Therefore, the wavefunctions in our case obey the relationship:

\begin{equation}
\psi_{R\boldsymbol{k}}(\boldsymbol{r}) = e^{-iR\boldsymbol{k}\cdot\boldsymbol{t}} \psi_{\boldsymbol{k}}(R^{-1}(\boldsymbol{r}-\boldsymbol{t}))
\end{equation}
or equivalently:
\begin{equation}
\psi_{R^{-1}\boldsymbol{k}+\boldsymbol{q}}(\boldsymbol{r}) = e^{i(R^{-1}\boldsymbol{k}+\boldsymbol{q})\cdot\boldsymbol{t}} \psi_{\boldsymbol{k}+R\boldsymbol{q}}(S\boldsymbol{r})
\end{equation}

Therefore, we must revise the symmetry relations as:
\f{eqnarray}{
g(\boldsymbol{k},S\boldsymbol{q})
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{S\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(S^{-1}\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(\boldsymbol{r}) \big> \\
&=& \big< \psi_{\boldsymbol{k}+S\boldsymbol{q}}(S\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{\boldsymbol{k}}(S\boldsymbol{r}) \big> \\
&=& e^{i(R^{-1}\boldsymbol{k}+\boldsymbol{q})\cdot\boldsymbol{t}} \big< \psi_{S^{-1}\boldsymbol{k}+\boldsymbol{q}}(\boldsymbol{r}) \big| \delta V_{\boldsymbol{q}}(\boldsymbol{r}) \big| \psi_{S^{-1}\boldsymbol{k}}(\boldsymbol{r}) \big> e^{-i(R^{-1}\boldsymbol{k})\cdot  \boldsymbol{t}}  \\
&=& e^{i\boldsymbol{q}\cdot \boldsymbol{t}} g(S^{-1}\boldsymbol{k},\boldsymbol{q})
\f}



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