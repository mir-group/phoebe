@page Theory Theory
@section WANHAR Wannier interpolation of electronic band structure

We assume that a third-party code has provided us with the single-particle hamiltonian in the wannier representation:

\begin{equation}
\langle 0n | H | R m \rangle
\end{equation}
where \f$R\f$ labels a bravais lattice vector in a supercell, \f$m\f$ and \f$n\f$ label two wannier functions.
A Wannier function is here denoted by the ket:
\begin{equation}
\langle R m \rangle
\end{equation}

The Wannier interpolation requires first to transform to reciprocal space:
\begin{equation}
H_{k,nm}^W = \sum_{R} e^{i k \cdot R} \langle 0n | H | R m \rangle
\end{equation}
This sum is typically speeded-up summing over irreducible lattice vectors:
\begin{equation}
H_{k,nm}^W = \sum_{R_{irr}} \frac{e^{i k \cdot R_{irr}} }{ d_{R_{irr}}} \langle 0n | H | R m \rangle
\end{equation}
where \f$d_{R_{irr}}\f$ is the degree of degeneracy of the irreducible bravais lattice vector \f$R_{irr}\f$.

This matrix is not diagonal, as we are working, typically, in the maximally localized wannier function gauge.
We thus diagonalize the matrix to pass to the Bloch representation:
\begin{equation}
H_{k,bb'}^B = [U_k^\dagger H_k^W U_k]_{bb'} = \delta_{bb'} \epsilon_{kb}
\end{equation}
Therefore finding the interpolated electronic energies at a generic \f$k\f$ point.

The velocity operator is computed with the Hellman-Feynman theorem, e.g. along the x direction as:
\begin{equation}
v^x_{k,bb'} = [U_k^\dagger \frac{H_(k+\delta_x)^W-H_(k-\delta_x)^W}{2 \delta_x} U_k]_{bb'}
\end{equation}
Whenever we find a set of degenerate bands at a given k point, we diagonalize the velocity operator in the degenerate subset, in order to uniquely define the velocity operator.



