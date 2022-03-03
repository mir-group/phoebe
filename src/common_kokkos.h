#ifndef COMMON_KOKKOS_H
#define COMMON_KOKKOS_H

#include <Kokkos_Core.hpp>

/** Define some useful classes for Kokkos-related calculations.
 */

// Kokkos View types for global memory spaces
using ComplexView1D = Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutRight>;
using ComplexView2D = Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutRight>;
using ComplexView3D = Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutRight>;
using ComplexView4D = Kokkos::View<Kokkos::complex<double> ****, Kokkos::LayoutRight>;
using ComplexView5D = Kokkos::View<Kokkos::complex<double> *****, Kokkos::LayoutRight>;
using IntView1D = Kokkos::View<int *, Kokkos::LayoutRight>;
using IntView2D = Kokkos::View<int **, Kokkos::LayoutRight>;
using DoubleView1D = Kokkos::View<double *, Kokkos::LayoutRight>;
using DoubleView2D = Kokkos::View<double **, Kokkos::LayoutRight>;
using DoubleView3D = Kokkos::View<double ***, Kokkos::LayoutRight>;
using DoubleView4D = Kokkos::View<double ****, Kokkos::LayoutRight>;
using DoubleView5D = Kokkos::View<double *****, Kokkos::LayoutRight>;

using HostComplexView1D = Kokkos::View<Kokkos::complex<double> *, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostComplexView2D = Kokkos::View<Kokkos::complex<double> **, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostComplexView3D = Kokkos::View<Kokkos::complex<double> ***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostComplexView4D = Kokkos::View<Kokkos::complex<double> ****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostComplexView5D = Kokkos::View<Kokkos::complex<double> *****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostDoubleView1D = Kokkos::View<double *, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostDoubleView2D = Kokkos::View<double **, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostDoubleView3D = Kokkos::View<double ***, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostDoubleView4D = Kokkos::View<double ****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;
using HostDoubleView5D = Kokkos::View<double *****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>;

// Kokkos Range types
using Range2D = Kokkos::MDRangePolicy<Kokkos::Rank<2, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;
using Range3D = Kokkos::MDRangePolicy<Kokkos::Rank<3, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;
using Range4D = Kokkos::MDRangePolicy<Kokkos::Rank<4, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;
using Range5D = Kokkos::MDRangePolicy<Kokkos::Rank<5, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;
using Range6D = Kokkos::MDRangePolicy<Kokkos::Rank<6, Kokkos::Iterate::Right, Kokkos::Iterate::Right>>;

/** This function does a batched diagonalization of M matrices A.
 * Each kernel works on a single matrix A
 *
 * @param A. On entry, a MxNxN tensor, identifying M hermitian matrices of size
 * NxN. matrix. On exit, A contains the eigenvectors of all the M matrices.
 * @param W: a MxN tensor, containing the N eigenvalues of each M-th matrix A.
 */
void kokkosZHEEV(ComplexView3D& A, DoubleView2D& W);

#endif
