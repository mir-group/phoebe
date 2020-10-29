@page Theory Theory
@section WANHAR Wannier interpolation of electronic band structure

A good review of the Wannier function formalism can be found [at this link] (https://arxiv.org/abs/1112.5411).

We assume that a third-party code has provided us with the single-particle hamiltonian in the wannier representation:

\begin{equation}
\langle \boldsymbol{0}n | H | \boldsymbol{R} m \rangle
\end{equation}
where \f$\boldsymbol{R}\f$ labels a bravais lattice vector in a supercell, \f$m\f$ and \f$n\f$ label two wannier functions.
A Wannier function is here denoted by the ket:
\begin{equation}
\langle \boldsymbol{R} m \rangle
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
v^x_{k,bb'} = [U_{\boldsymbol{k}}^\dagger \frac{H_{\boldsymbol{k}+\boldsymbol{\delta}_x}^W-H_{\boldsymbol{k}-\boldsymbol{\delta}_x}^W}{2 \delta_x} U_{\boldsymbol{k}}]_{bb'}
\end{equation}
Whenever we find a set of degenerate bands at a given k point, we diagonalize the velocity operator in the degenerate subset, in order to uniquely define the velocity operator.



