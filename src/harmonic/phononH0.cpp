#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <math.h>
#include <iostream>
#include <complex>

#include "constants.h"
#include "phononH0.h"
#include "exceptions.h"

Eigen::MatrixXd PhononH0::wsinit(const Eigen::Matrix3d& unitCell) {
	const int nx=2;
	int index = 0;
	const int nrwsx = 200;

	Eigen::MatrixXd tmpResult(3,nrwsx);

	for ( int ir=-nx; ir<=nx; ir++ ) {
		for ( int jr=-nx; jr<=nx; jr++ ) {
			for ( int kr=-nx; kr<=nx; kr++ ) {
				for ( int i=0; i<3; i++ ) {
					tmpResult(i,index) = unitCell(i,0)*ir
							+ unitCell(i,1)*jr + unitCell(i,2)*kr;
				}

				if ( tmpResult.col(index).transpose()*tmpResult.col(index) > 1.0e-6 ) {
					index += 1;
				}
				if ( index > nrwsx ) {
					error("WSInit > nrwsx",1);
				}
			}
		}
	}
	int nrws = index;

	Eigen::MatrixXd result(3,nrws);
	for ( int i=0; i<nrws; i++ ) {
		result.col(i) = tmpResult.col(i);
	}
	return result;
}

double PhononH0::wsweight(const Eigen::VectorXd& r,
		const Eigen::MatrixXd& rws) {
	//! wsweights assigns this weight:
	//! - if a point is inside the Wigner-Seitz cell:    weight=1
	//! - if a point is outside the WS cell:             weight=0
	//! - if a point q is on the border of the WS cell, it finds the number N
	//!   of translationally equivalent point q+G  (where G is a lattice vector)
	//!   that are also on the border of the cell. Then: weight = 1/N
	//
	//! I.e. if a point is on the surface of the WS cell of a cubic lattice
	//! it will have weight 1/2; on the vertex of the WS it would be 1/8;
	//! the K point of an hexagonal lattice has weight 1/3 and so on.

	// rws: contains the list of nearest neighbor atoms
	// r: the position of the reference point
	// rws.cols(): number of nearest neighbors

	int nreq = 1;
	double rrt, ck;

	for ( int ir=0; ir<rws.cols(); ir++ ) {
		rrt = r.transpose() * rws.col(ir);
		ck = rrt - (rws.col(ir).transpose() * rws.col(ir)).value() / 2.;
		if ( ck > 1.0e-6 ) {
			return 0.;
		}
		if ( abs(ck) < 1.0e-6 ) {
			nreq += 1;
		}
	}

	double x = 1. / (double)nreq;

	return x;
}

void PhononH0::longRangeTerm(Eigen::Tensor<std::complex<double>,4>& dyn,
		const Eigen::VectorXd& q,
		const int sign) {
	// this subroutine is the analogous of rgd_blk in QE
	// compute the rigid-ion (long-range) term for q
	// The long-range term used here, to be added to or subtracted from the
	// dynamical matrices, is exactly the same of the formula introduced in:
	// X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is
	// implemented: the Ewald parameter alpha must be large enough to
	// have negligible r-space contribution

	//   complex(DP) :: dyn(3,3,numAtoms,numAtoms) ! dynamical matrix
	//   real(DP) &
	//        q(3),           &! q-vector
	//        sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term

	// alph is the Ewald parameter, geg is an estimate of G^2
	// such that the G-space sum is convergent for that alph
	// very rough estimate: geg/4/alph > gmax = 14
	// (exp (-14) = 10^-6)

	double geg, gp2, r; //  <q+G| dielectricMatrix | q+G>,  For 2d loto: gp2, r
	int nr1x, nr2x, nr3x;
	Eigen::VectorXd zag(3), zbg(3), zcg(3), fnat(3);
	Eigen::MatrixXd reff(2,2);
	double fac, facgd, arg;
	std::complex<double> facg;

	double gmax = 14.;
	double alph = 1.;
	geg = gmax * alph * 4.;

	// Estimate of nr1x,nr2x,nr3x generating all vectors up to G^2 < geg
	// Only for dimensions where periodicity is present, e.g. if nr1=1
	// and nr2=1, then the G-vectors run along nr3 only.
	// (useful if system is in vacuum, e.g. 1D or 2D)

	if ( qCoarseGrid(0) == 1 ) {
		nr1x = 0;
	} else {
		nr1x = (int) ( sqrt(geg) /
				sqrt( reciprocalUnitCell.col(0).transpose() * reciprocalUnitCell.col(0) )) + 1;
	}
	if ( qCoarseGrid(1) == 1 ) {
		nr2x = 0;
	} else {
		nr2x = (int) ( sqrt(geg) /
				sqrt( reciprocalUnitCell.col(1).transpose() * reciprocalUnitCell.col(1) )) + 1;
	}
	if ( qCoarseGrid(2) == 1 ) {
		nr3x = 0;
	} else {
		nr3x = (int) ( sqrt(geg) /
				sqrt( reciprocalUnitCell.col(2).transpose() * reciprocalUnitCell.col(2) )) + 1;
	}

	if ( abs(sign) != 1. ) {
		error("wrong value for sign", 1);
	}

	if ( loto_2d ) {
		fac = sign * e2 * fourPi / volumeUnitCell * 0.5 * latticeParameter / reciprocalUnitCell(3,3);
		reff.setZero();
		for ( int i=0; i<2; i++ ) {
			for ( int j=0; j<2; j++ ) {
				reff(i,j) = dielectricMatrix(i,j) * 0.5 * twoPi / reciprocalUnitCell(3,3); // (eps)*c/2 in 2pi/a units
			}
		}
		for ( int i=0; i<2; i++ ) {
			reff(i,i) = reff(i,i) - 0.5 * twoPi / reciprocalUnitCell(3,3); // (-1)*c/2 in 2pi/a units
		}
	} else {
		fac = sign * e2 * fourPi / volumeUnitCell;
	}

	Eigen::VectorXd g(3);
	std::complex<double> phase;
	for ( int m1=-nr1x; m1<=nr1x; m1++ ) {
		for ( int m2=-nr2x; m2<=nr2x; m2++ ) {
			for ( int m3=-nr3x; m3<=nr3x; m3++ ) {
				g(0) = m1*reciprocalUnitCell(0,0) + m2*reciprocalUnitCell(0,1) + m3*reciprocalUnitCell(0,2);
				g(1) = m1*reciprocalUnitCell(1,0) + m2*reciprocalUnitCell(1,1) + m3*reciprocalUnitCell(1,2);
				g(2) = m1*reciprocalUnitCell(2,0) + m2*reciprocalUnitCell(2,1) + m3*reciprocalUnitCell(2,2);

				if (loto_2d) {
					geg = (g.transpose()*g).value();
					r = 0.;
					gp2 = g(0)*g(0) + g(1)*g(1);
					if ( gp2 > 1.0e-8 ) {
						r = g(0) * reff(0,0) * g(0) + g(0) * reff(0,1) * g(1) + g(1) * reff(1,0) * g(0) + g(1) * reff(1,1) * g(1);
						r = r / gp2;
					}
				} else {
					geg = (g.transpose() * dielectricMatrix * g).value();
				}

				if ( geg > 0. && geg / alph / 4. < gmax ) {

					if ( loto_2d ) {
						facgd = fac * exp( - geg / alph / 4.) / sqrt(geg) / ( 1. + r * sqrt(geg) );
					} else {
						facgd = fac * exp( - geg / alph / 4. ) / geg;
					}

					for ( int na=0; na<numAtoms; na++ ) {

						for ( int i=0; i<3; i++ ) {
							zag(i) = g(0) * bornCharges(na,0,i) + g(1) * bornCharges(na,1,i) + g(2) * bornCharges(na,2,i);
							fnat(i) = 0.;
							for ( int nb=0; nb<numAtoms; nb++ ) {
								arg = ( (atomicPositions.row(na)-atomicPositions.row(nb)) * g ).value();
								arg *= twoPi;
								zcg(i) = g(0) * bornCharges(nb,0,i) + g(1) * bornCharges(nb,1,i) + g(2) * bornCharges(nb,2,i);
								fnat(i) += zcg(i) * cos(arg);
							}
						}

						for ( int i=0; i<3; i++ ) {
							for ( int j=0; j<3; j++ ) {
								dyn(i,j,na,na) += - facgd * zag(i) * fnat(j);
							}
						}
					}
				}

				g += q;

				if ( loto_2d ) {
					geg = (g.transpose()*g).value();
					r = 0.;
					gp2 = g(0)*g(0) + g(1)*g(1);
					if ( gp2 > 1.0e-8) {
						r = g(0) * reff(0,0)*g(0) + g(0) * reff(0,1)*g(1) + g(1) * reff(1,0)*g(0) + g(1) * reff(1,1)*g(1);
						r = r / gp2;
					}
				} else {
					geg = (g.transpose() * dielectricMatrix * g).value();
				}

				if ( geg > 0. && geg / alph / 4. < gmax ) {

					if ( loto_2d ) {
						facgd = fac * exp(-geg/alph/4.) / sqrt(geg) / (1.+r*sqrt(geg));
					} else {
						facgd = fac * exp( - geg / alph / 4. ) / geg;
					}

					for ( int nb=0; nb<numAtoms; nb++ ) {
						for ( int i=0; i<3; i++ ) {
							zbg(i) = g(0) * bornCharges(nb,0,i) + g(1) * bornCharges(nb,1,i) + g(2) * bornCharges(nb,2,i);
						}
						for ( int na=0; na<numAtoms; na++ ) {
							for ( int i=0; i<3; i++ ) {
								zag(i) = g(0) * bornCharges(na,0,i) + g(1) * bornCharges(na,1,i) + g(2) * bornCharges(na,2,i);
							}
							arg = ( (atomicPositions.row(na)-atomicPositions.row(nb)) * g ).value();
							arg *= twoPi;
							phase = {cos(arg), sin(arg)};
							facg = facgd * phase;
							for ( int i=0; i<3; i++ ) {
								for ( int j=0; j<3; j++ ) {
									dyn(i,j,na,nb) += facg * zag(i) * zbg(j);
								}
							}
						}
					}
				}
			}
		}
	}
}

void PhononH0::nonAnaliticTerm(const Eigen::VectorXd& q,
		Eigen::Tensor<std::complex<double>,4>& dyn ) {
	// add the nonAnalytical term with macroscopic electric fields
	// numAtoms: number of atoms in the cell (in the supercell in the case
	//       of a dyn.mat. constructed in the mass approximation)
	// atomicSpecies(na): atom in the original cell corresponding to
	//                atom na in the supercell
	// dyn(3,3,numAtoms,numAtoms) ! dynamical matrix
	// q(3), & ! polarization vector
	// dielectricMatrix(3,3),              & ! dielectric constant tensor
	// bornCharges(3,3,numAtoms),        & ! effective charges tensor
	// volumeUnitCell                      ! unit cell volume

	double qeq = (q.transpose() * dielectricMatrix * q).value();
	if ( qeq < 1.e-8 ) {
		warning("A direction for q was not specified: "
				"TO-LO splitting will be absent");
		return;
	}

	Eigen::VectorXd zag(3), zbg(3);

	int iType, jType;

	for ( int it=0; it<numAtoms; it++ ) {
		iType = atomicSpecies(it);
		for ( int jt=0; jt<numAtoms; jt++ ) {
			jType = atomicSpecies(jt);

			for ( int i=0; i<3; i++ ) {
				zag(i) = q(0)*bornCharges(iType,0,i) +  q(1)*bornCharges(iType,1,i) +
						q(2)*bornCharges(iType,2,i);
				zbg(i) = q(0)*bornCharges(jType,0,i) +  q(1)*bornCharges(jType,1,i) +
						q(2)*bornCharges(jType,2,i);
			}

			for ( int i=0; i<3; i++ ) {
				for ( int j=0; j<3; j++ ) {
					dyn(i,j,it,jt) += fourPi * e2 * zag(i) * zbg(j) / qeq / volumeUnitCell;
				}
			}
		}
	}
}

void PhononH0::nonAnalIFC(const Eigen::VectorXd& q,
		Eigen::Tensor<std::complex<double>, 4>& f_of_q) {
	//     add the nonanalytical term with macroscopic electric fields

	// !  numAtoms: number of atoms in the cell (in the supercell in the case
	// !       of a dyn.mat. constructed in the mass approximation)
	// !  atomicSpecies(na): atom in the original cell corresponding to
	// !                atom na in the supercell

	//	complex(DP), intent(inout) :: dyn(3,3,numAtoms,numAtoms),f_of_q(3,3,numAtoms,numAtoms) ! dynamical matrix
	// real(DP), intent(in) :: q(3),  &! polarization vector
	//      &       dielectricMatrix(3,3),     &! dielectric constant tensor
	//      &       bornCharges(3,3,numAtoms),   &! effective charges tensor
	//      &       volumeUnitCell            ! unit cell volume

	Eigen::VectorXd zag(3), zbg(3); // eff. charges  times g-vector

	if ( q.transpose() * q == 0. ) {
		return;
	}

	double qeq = q.transpose() * dielectricMatrix * q;

	if ( qeq < 1.0e-8 ) {
		return;
	}

	int na_blk, nb_blk;

	double factor = fourPi * e2 / qeq / volumeUnitCell
			/ (qCoarseGrid(0)*qCoarseGrid(1)*qCoarseGrid(2));

	for ( int na=0; na<numAtoms; na++ ) {
		na_blk = atomicSpecies(na);
		for ( int nb=0; nb<numAtoms; nb++ ) {
			nb_blk = atomicSpecies(nb);
			for ( int i=0; i<3; i++ ) {
				zag(i) = q(0)*bornCharges(na_blk,0,i) +  q(1)*bornCharges(na_blk,1,i) +
						q(2)*bornCharges(na_blk,2,i);
				zbg(i) = q(0)*bornCharges(nb_blk,0,i) +  q(1)*bornCharges(nb_blk,1,i) +
						q(2)*bornCharges(nb_blk,2,i);
			}
			for ( int i=0; i<3; i++ ) {
				for ( int j=0; j<3; j++ ) {
					f_of_q(i,j,na,nb) = factor * zag(i) * zbg(j);
				}
			}
		}
	}
}

void PhononH0::shortRangeTerm(Eigen::Tensor<std::complex<double>, 4>& dyn,
		const Eigen::VectorXd& q,
		Eigen::Tensor<std::complex<double>, 4>& f_of_q) {
	// calculates the dynamical matrix at q from the (short-range part of the)
	// force constants21

	Eigen::VectorXd r(3), r_ws(3);
	double total_weight, arg, weight;
	static bool first = true;

	int n1ForCache, n2ForCache, n3ForCache;

	const int nr1Big = 2 * qCoarseGrid(0);
	const int nr2Big = 2 * qCoarseGrid(1);
	const int nr3Big = 2 * qCoarseGrid(2);

	static Eigen::Tensor<double, 5> wscache(2*nr3Big+1, 2*nr2Big+1, 2*nr1Big+1,
			numAtoms, numAtoms);

	double x;
	if ( first ) {
		first = false;
		for ( int na=0; na<numAtoms; na++ ) {
			for ( int nb=0; nb<numAtoms; nb++ ) {
				total_weight = 0.;

				// sum over r vectors in the supercell - very safe range!

				for ( int n1=-nr1Big; n1<=nr1Big; n1++ ) {
					n1ForCache = n1 + nr1Big;
					for ( int n2=-nr2Big; n2<=nr2Big; n2++ ) {
						n2ForCache = n2 + nr2Big;
						for ( int n3=-nr3Big; n3<=nr3Big; n3++ ) {
							n3ForCache = n3 + nr3Big;

							for ( int i=0; i<3; i++ ) {
								r(i) = n1 * directUnitCell(i,0) + n2 * directUnitCell(i,1) + n3*directUnitCell(i,2);
								r_ws(i) = r(i) + atomicPositions(na,i) - atomicPositions(nb,i);
								if ( frozenPhonon ) {
									r_ws(i) = r(i) + atomicPositions(nb,i) - atomicPositions(na,i);
								}
							}

							x = wsweight(r_ws, rws);

							wscache(n3ForCache, n2ForCache, n1ForCache, nb, na)
							= x;
							total_weight += x;
						}
					}
				}

				if ( abs( total_weight - qCoarseGrid(0) *
						qCoarseGrid(1) * qCoarseGrid(2) ) > 1.0e-8 ) {
					error("wrong total_weight", 1);
				}
			}
		}
	} // first_time only

	int m1, m2, m3;
	std::complex<double> phase;

	for ( int na=0; na<numAtoms; na++ ) {
		for ( int nb=0; nb<numAtoms; nb++ ) {
			for ( int n1=-nr1Big; n1<nr1Big; n1++ ) {
				n1ForCache = n1 + nr1Big;
				for ( int n2=-nr2Big; n2<nr2Big; n2++ ) {
					n2ForCache = n2 + nr2Big;
					for ( int n3=-nr3Big; n3<nr3Big; n3++ ) {
						n3ForCache = n3 + nr3Big;

						// sum over r vectors in the supercell, very safe range
						for ( int i=0; i<3; i++ ) {
							r(i) = n1 * directUnitCell(i,0)
							+ n2 * directUnitCell(i,1)
							+ n3 * directUnitCell(i,2);
						}

						weight = wscache(n3ForCache,n2ForCache,n1ForCache,nb,na);
						if ( weight > 0. ) {

							// find vector corresponding to r in original cell

							m1 = ( n1 + 1 ) % qCoarseGrid(0);
							if ( m1 <= 0 ) { m1 += qCoarseGrid(0); };
							m2 = ( n2 + 1 ) % qCoarseGrid(1);
							if ( m2 <= 0 ) { m2 += qCoarseGrid(1); };
							m3 = ( n3 + 1 ) % qCoarseGrid(2);
							if ( m3 <= 0 ) { m3 += qCoarseGrid(2); };
							m1 += -1;
							m2 += -1;
							m3 += -1;

							// FOURIER TRANSFORM

							// note: maybe send m1 in m1-1 and m2,m3

							arg = twoPi * (q.transpose() * r).value();
							phase = {cos(arg),-sin(arg)};

							for ( int ipol=0; ipol<3; ipol++ ) {
								for ( int jpol=0; jpol<3; jpol++ ) {
									dyn(ipol,jpol,na,nb) +=
											(forceConstants(m1,m2,m3,ipol,jpol,na,nb)
													+ f_of_q(ipol,jpol,na,nb))
													* phase * weight;
								}
							}

						}
					}
				}
			}
		}
	}
}

void PhononH0::dyndiag(Eigen::Tensor<std::complex<double>,4>& dyn,
		Eigen::VectorXd& energies,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors) {
	// diagonalise the dynamical matrix
	// On input:  speciesMasses = masses, in amu
	// On output: w2 = energies, z = displacements

	// fill the two-indices dynamical matrix

	int iType, jType;
	std::complex<double> cx;
	Eigen::MatrixXcd dyn2Tmp(numBands, numBands);
	Eigen::MatrixXcd dyn2(numBands, numBands);

	for (int iat = 0; iat<numAtoms; iat++) {
		for (int jat = 0; jat<numAtoms; jat++) {
			for (int ipol = 0; ipol<3; ipol++) {
				for (int jpol = 0; jpol<3; jpol++) {
					dyn2Tmp(iat*3 +ipol, jat*3 +jpol) = dyn(ipol,jpol,iat,jat);
				}
			}
		}
	}

	// impose hermiticity

	dyn2 = dyn2Tmp + dyn2Tmp.adjoint();
	dyn2 *= 0.5;

	//  divide by the square root of masses

	for ( int iat=0; iat<numAtoms; iat++ ) {
		iType = atomicSpecies(iat);
		for ( int jat=0; jat<numAtoms; jat++ ) {
			jType = atomicSpecies(jat);
			for ( int ipol=0; ipol<3; ipol++ ) {
				for ( int jpol=0; jpol<3; jpol++ ) {
					dyn2(iat*3 + ipol, jat*3 + jpol) /=
							sqrt(speciesMasses(iType)*speciesMasses(jType));
				}
			}
		}
	}

	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> eigensolver(dyn2);

	Eigen::VectorXd w2 = eigensolver.eigenvalues();

	for ( int i=0; i<numBands; i++ ) {
		if ( energies(i) < 0 ) {
			energies(i) = sqrt(-w2(i));
		} else {
			energies(i) = sqrt(w2(i));
		}
	}
	std::cout << energies.transpose() << "\n";
	std::cout << energies.transpose() * ryToCmm1 << "\n";

	Eigen::MatrixXcd zTemp = eigensolver.eigenvectors();

	//  displacements are eigenvectors divided by sqrt(speciesMasses)

	for ( int iband=0; iband<numBands; iband++ ) {
		for ( int iat=0; iat<numAtoms; iat++ ) {
			iType = atomicSpecies(iat);
			for ( int ipol=0; ipol<3; ipol++ ) {
				eigenvectors(ipol,iat,iband) = zTemp(iat*3 + ipol, iband)
								/ sqrt(speciesMasses(iType));
			}
		}
	}
};

// TODO: we will move this functionality somewhere else with a kpoints class
void cryst_to_cart(Eigen::VectorXd& vec, const Eigen::MatrixXd& trmat,
		const int iflag) {
	//  !     This routine transforms the atomic positions or the k-point
	//  !     components from crystallographic to cartesian coordinates
	//  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
	//  !     Output cartesian coordinates are stored in the input ('vec') array
	//  integer, intent(in) :: nvec, iflag
	//  ! nvec:  number of vectors (atomic positions or k-points)
	//  !        to be transformed from crystal to cartesian and vice versa
	//  ! iflag: gives the direction of the transformation
	//  real(DP), intent(in) :: trmat (3, 3)
	//  ! trmat: transformation matrix
	//  ! if iflag=1:
	//  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
	//  !          = bg ,  basis of the reciprocal-space lattice, for k-points
	//  ! if iflag=-1: the opposite
	//  real(DP), intent(inout) :: vec (3, nvec)
	//  ! coordinates of the vector (atomic positions or k-points) to be
	//  ! transformed - overwritten on output
	//  !
	//  !    local variables
	//  !
	//  integer :: nv, kpol
	//  ! counter on vectors
	//  ! counter on polarizations
	//  real(DP) :: vau (3)
	//  ! workspace
	//  !
	//  !     Compute the cartesian coordinates of each vectors
	//  !     (atomic positions or k-points components)
	//  !
	if ( iflag == 1 ) {
		for ( int kpol=0; kpol<3; kpol++ )
		{
			vec = trmat * vec;
		}
	} else {
		for ( int kpol=0; kpol<3; kpol++ )
		{
			vec = vec * trmat;
		}
	}
}


void PhononH0::setup(
		const Eigen::MatrixXd& directUnitCell_,
		const Eigen::MatrixXd& reciprocalUnitCell_,
		const double& latticeParameter_,
		const double& volumeUnitCell_,
		const Eigen::MatrixXi& atomicSpecies_,
		const Eigen::VectorXd& speciesMasses_,
		const Eigen::MatrixXd& atomicPositions_,
		const Eigen::MatrixXd& dielectricMatrix_,
		const Eigen::Tensor<double, 3>& bornCharges_,
		Eigen::VectorXi& qCoarseGrid_,
		const Eigen::Tensor<double, 7> forceConstants_) {

	// in this section, we save as class properties a few variables
	// that are needed for the diagonalization of phonon frequencies

	// this variables might be used for extending future functionalities.
	// for the first tests, they can be left at these default values
	// in the future, we might expose them to the user input
	na_ifc = false;
	loto_2d = false;
	frozenPhonon = false;

	if ( dielectricMatrix.sum() > 0. ) {
		hasDielectric = true;
	} else {
		hasDielectric = false;
	}

	numAtoms = atomicPositions.rows();
	numBands = numAtoms * 3;

	directUnitCell = directUnitCell_;
	reciprocalUnitCell = reciprocalUnitCell_;
	latticeParameter = latticeParameter_;
	volumeUnitCell = volumeUnitCell_;
	atomicSpecies = atomicSpecies_;
	speciesMasses = speciesMasses_;
	atomicPositions = atomicPositions_;
	dielectricMatrix = dielectricMatrix_;
	bornCharges = bornCharges_;
	qCoarseGrid_ = qCoarseGrid_;
	forceConstants = forceConstants_;

	// now, I initialize an auxiliary set of vectors that are needed
	// for the diagonalization, which are precomputed once and for all.

	Eigen::Matrix3d directUnitCellSup(3,3);
	directUnitCellSup.col(0) = directUnitCell.col(0) * qCoarseGrid(0);
	directUnitCellSup.col(1) = directUnitCell.col(1) * qCoarseGrid(1);
	directUnitCellSup.col(2) = directUnitCell.col(2) * qCoarseGrid(2);
	rws = PhononH0::wsinit(directUnitCellSup);
}

void PhononH0::diagonalize(const Eigen::VectorXd& q,
		Eigen::VectorXd& energies,
		Eigen::Tensor<std::complex<double>,3>& eigenvectors) {
	// to be executed at every qpoint to get phonon frequencies and wavevectors

	Eigen::Tensor<std::complex<double>, 4> dyn(3,3,numAtoms,numAtoms);
	dyn.setZero();

	Eigen::Tensor<std::complex<double>,4> f_of_q(3,3,numAtoms,numAtoms);
	f_of_q.setZero();

	Eigen::VectorXd qhat(3);
	double qq;

	// for now, this part is not executed
	if ( na_ifc ) {
		qq = sqrt( q.transpose()*q ); // q is the qpoint coordinate
		if ( abs(qq) < 1.0e-8 ) {
			qq = 1.;
		}

		qhat = q / qq;

		nonAnalIFC(qhat, f_of_q);
	}

	// first, the short range term, which is just a Fourier transform
	shortRangeTerm(dyn, q, f_of_q);

	// then the long range term, which uses some convergence
	// tricks by X. Gonze et al.
	if ( hasDielectric && !na_ifc ) {
		longRangeTerm(dyn, q, +1.);
	}

	// finally, the nonanalytic term from Born charges
	if ( !loto_2d && na_ifc ) {
		qhat = q.transpose() * directUnitCell;
		if ( abs( qhat(0) - round(qhat(0) ) ) <= 1.0e-6 &&
				abs( qhat(1) - round(qhat(1) ) ) <= 1.0e-6 &&
				abs( qhat(2) - round(qhat(2) ) ) <= 1.0e-6 ) {
			// q = 0 : we need the direction q => 0 for the non-analytic part

			qq = sqrt( ( qhat.transpose()*qhat ).value() );
			if (qq != 0. ) {
				qhat /= qq;
			}
			nonAnaliticTerm(qhat, dyn);
		}
	}

	// once everything is ready, here we scale by masses and diagonalize
	dyndiag(dyn, energies, eigenvectors);
};


//void PhononH0::setAcousticSumRule(const std::string sumRule) {
//	double norm2;
//	//INTEGER :: axis, n, i, j, na, nb, n1,n2,n3, m,p,k,l,q,r, i1,j1,na1
//	//  type vector
//	//     real(DP),pointer :: vec(:,:,:,:,:,:,:)
//	//  end type vector
//	//  type (vector) u(6*3*nat)
//	//  ! These are the "vectors" associated with the sum rules on force-constants
//
//	//  integer :: u_less(6*3*nat),n_less,i_less
//	//  ! indices of the vectors u that are not independent to the preceding ones,
//	//  ! n_less = number of such vectors, i_less = temporary parameter
//	//  !
//	//  integer, allocatable :: ind_v(:,:,:)
//	//  real(DP), allocatable :: v(:,:)
//	//  ! These are the "vectors" associated with symmetry conditions, coded by
//	//  ! indicating the positions (i.e. the seven indices) of the non-zero elements (there
//	//  ! should be only 2 of them) and the value of that element. We do so in order
//	//  ! to limit the amount of memory used.
//	//  !
//	//  real(DP), allocatable :: w(:,:,:,:,:,:,:), x(:,:,:,:,:,:,:)
//	//  ! temporary vectors and parameters
//	//  real(DP) :: scal,norm2, sum
//	//  !
//	//  real(DP) :: zeu_u(6*3,3,3,nat)
//	//  ! These are the "vectors" associated with the sum rules on effective charges
//	//  !
//	//  integer :: zeu_less(6*3),nzeu_less,izeu_less
//	//  ! indices of the vectors zeu_u that are not independent to the preceding ones,
//	//  ! nzeu_less = number of such vectors, izeu_less = temporary parameter
//	//  !
//	//  real(DP) :: zeu_w(3,3,nat), zeu_x(3,3,nat)
//	//  ! temporary vectors
//
//	// Initialization. n is the number of sum rules to be considered
//	// (if sumRule!="simple")
//	// and 'axis' is the rotation axis in the case of a 1D system (i.e. the
//	// rotation axis is (Ox) if axis='1', (Oy) if axis='2' and (Oz) if axis='3')
//
//	if ( ( sumRule != "simple" ) && ( sumRule != "crystal" ) ) {
//		error("invalid Acoustic Sum Rule", 1);
//	}
//	//  if ( ( sumRule != "simple" ) && ( sumRule != "crystal" ) &&
//	//	   ( sumRule != "one-dim") && ( sumRule != "zero-dim" ) ) {
//	//     error("invalid Acoustic Sum Rule", 1);
//	//  }
//
//	if ( sumRule == "simple" ) {
//
//		// Simple Acoustic Sum Rule on effective charges
//
//		double sum;
//
//		for ( int i=0; i<3; i++ ) {
//			for ( int j=0; j<3; j++ ) {
//				sum = 0.;
//				for ( int na=0; na<numAtoms; na++ ) {
//					sum += bornCharges(i,j,na);
//				}
//				for ( int na=0; na<3; na++ ) {
//					bornCharges(na,i,j) -= sum / numAtoms;
//				}
//			}
//		}
//
//		// Simple Acoustic Sum Rule on force constants in real space
//
//		for ( int i=0; i<3; i++ ) {
//			for ( int j=0; j<3; j++ ) {
//				for ( int na=0; na<numAtoms; na++ ) {
//					sum = 0.;
//					for ( int nb=0; nb<numAtoms; nb++ ) {
//
//						for ( int n1=0; n1<qCoarseGrid(0); n1++ ) {
//							for ( int n2=0; n2<qCoarseGrid(1); n2++ ) {
//								for ( int n3=0; n3<qCoarseGrid(2); n3++ ) {
//									sum += forceConstants(n1,n2,n3,i,j,na,nb);
//								}
//							}
//						}
//					}
//					forceConstants(1,1,1,i,j,na,na) -= sum;
//				}
//			}
//		}
//	} else {
//
//		//  if ( sumRule == "crystal") n=3
//		//  if ( sumRule == "one-dim" ) {
//		//     // the direction of periodicity is the rotation axis
//		//     // It will work only if the crystal axis considered is one of
//		//     // the cartesian axis (typically, ibrav=1, 6 or 8, or 4 along the
//		//     // z-direction)
//		//     if ( qCoarseGrid(0)*qCoarseGrid(1)*qCoarseGrid(2) == 1 ) axis = 3;
//		//     if ( (qCoarseGrid(0)!=1) && (qCoarseGrid(1)*qCoarseGrid(2)==1)) axis = 1;
//		//     if ( (qCoarseGrid(1)!=1) && (qCoarseGrid(0)*qCoarseGrid(2)==1)) axis = 2;
//		//     if ( (qCoarseGrid(2)!=1) && (qCoarseGrid(0)*qCoarseGrid(1)==1)) axis = 3;
//		//     if ( ((qCoarseGrid(0)!=1) && (qCoarseGrid(1)!=1))
//		//    		 || ((qCoarseGrid(1)!=1) && (qCoarseGrid(2)!=1))
//		//    		 || ((qCoarseGrid(0)!=1) && (qCoarseGrid(2) != 1) )) {
//		//        error("too many directions of periodicity in 1D system", 2);
//		//     }
//		//     if ( (ibrav.ne.1).and.(ibrav.ne.6).and.(ibrav.ne.8).and. &
//		//          ((ibrav.ne.4).or.(axis.ne.3)) ) then
//		//        write(stdout,*) 'sumRule: rotational axis may be wrong'
//		//     endif
//		//     write(stdout,'("sumRule rotation axis in 1D system= ",I4)') axis
//		//     n=4
//		//  }
//
//		//  if ( sumRule == "zero-dim") n=6;
//
//
//
//		// Acoustic Sum Rule on effective charges
//
//		// generating the vectors of the orthogonal of the subspace to project
//		// the effective charges matrix on
//
//		Eigen::Tensor<double,4> zeu_u(6*3,3,3,numAtoms);
//		zeu_u.setZero();
//		Eigen::Tensor<double,3> zeu_new = bornCharges;
//
//		int r;
//		int p = 0;
//		for ( int i=0; i<3; i++ ) {
//			for ( int j=0; j<3; j++ ) {
//				for ( int iat=0; iat<numAtoms; iat++ ) {
//					// These are the 3*3 vectors associated with the
//					// translational acoustic sum rules
//					zeu_u(p,i,j,iat) = 1.;
//				}
//				p += 1;
//			}
//		}
//
//		//  if (n.eq.4) then
//		//     do i=1,3
//		//        ! These are the 3 vectors associated with the
//		//        ! single rotational sum rule (1D system)
//		//        p=p+1
//		//        do na=1,nat
//		//           zeu_u(p,i,MOD(axis,3)+1,na)=-tau(MOD(axis+1,3)+1,na)
//		//           zeu_u(p,i,MOD(axis+1,3)+1,na)=tau(MOD(axis,3)+1,na)
//		//        enddo
//		//        !
//		//     enddo
//		//  endif
//		//  !
//		//  if (n.eq.6) then
//		//     do i=1,3
//		//        do j=1,3
//		//           ! These are the 3*3 vectors associated with the
//		//           ! three rotational sum rules (0D system - typ. molecule)
//		//           p=p+1
//		//           do na=1,nat
//		//              zeu_u(p,i,MOD(j,3)+1,na)=-tau(MOD(j+1,3)+1,na)
//		//              zeu_u(p,i,MOD(j+1,3)+1,na)=tau(MOD(j,3)+1,na)
//		//           enddo
//		//           !
//		//        enddo
//		//     enddo
//		//  endif
//
//		// Gram-Schmidt orthonormalization of the set of vectors created.
//
//		// temporary vectors
//		Eigen::Tensor<double,3> zeu_w(3,3,numAtoms), zeu_x(3,3,numAtoms);
//		Eigen::Tensor<double,3> tempZeu(3,3,numAtoms);
//		zeu_w.setZero();
//		zeu_x.setZero();
//		Eigen::VectorXi zeu_less(6*3);
//
//		double scal;
//
//		Eigen::array<long,4> offset4;
//		Eigen::array<long,4> extent4;
//
//		int nzeu_less = 0;
//		for ( int k=0; k<p; k++ ) {
//			offset4 = {k,0,0,0};
//			extent4 = {k+1,3,3,numAtoms};
//			zeu_w = zeu_u.slice(offset4, extent4);
//			zeu_x = zeu_u.slice(offset4, extent4);
//			// equivalent to:
//			//zeu_w(:,:,:) = zeu_u(k,:,:,:);
//			//zeu_x(:,:,:) = zeu_u(k,:,:,:);
//			for ( int q=0; q<k-1; q++ ) {
//				r = 1;
//				for ( int izeu_less=0; izeu_less<nzeu_less; izeu_less++ ) {
//					if ( zeu_less(izeu_less) == q ) { r = 0; };
//				}
//				if ( r != 0 ) {
//					offset4 = {q,0,0,0};
//					extent4 = {q+1,3,3,numAtoms};
//					tempZeu = zeu_u.slice(offset4,extent4);
//					// i.e. zeu_u(q,:,:,:)
//					sp_zeu(zeu_x, tempZeu, scal);
//					zeu_w -= scal * tempZeu;
//				}
//			}
//			sp_zeu(zeu_w, zeu_w, norm2);
//			if ( norm2 > 1.0e-16 ) {
//				offset4 = {k,0,0,0};
//				extent4 = {k+1,3,3,numAtoms};
//				// zeu_u(k,:,:,:)
//				zeu_u.slice(offset4,extent4) = zeu_w / sqrt(norm2);
//			} else {
//				zeu_less(nzeu_less) = k;
//				nzeu_less += 1;
//			}
//		}
//
//		// Projection of the effective charge "vector" on the orthogonal of the
//		// subspace of the vectors verifying the sum rules
//
//		zeu_w.setZero();
//		for ( int k=0; k<p; k++ ) {
//			r = 1;
//			for ( int izeu_less=0; izeu_less<nzeu_less; izeu_less++ ) {
//				if ( zeu_less(izeu_less) == k ) { r = 0; };
//			}
//			if ( r != 0 ) {
//				// zeu_u(k,:,:,:)
//				offset4 = {k,0,0,0};
//				extent4 = {k+1,3,3,numAtoms};
//				// zeu_u(k,:,:,:)
//				zeu_x = zeu_u.slice(offset4,extent4);
//				sp_zeu(zeu_x, zeu_new, scal);
//				zeu_w += scal * zeu_u.slice(offset4,extent4);
//			}
//		}
//
//		// Final substraction of the former projection to the initial zeu, to get
//		// the new "projected" zeu
//
//		zeu_new -= zeu_w;
//		sp_zeu(zeu_w, zeu_w, norm2);
//		std::cout << "Norm of the difference between old and new effective charges: "
//				<< sqrt(norm2);
//		bornCharges = zeu_new;
//
//		// Acoustic Sum Rule on force constants
//
//		// generating the vectors of the orthogonal of the subspace to project
//		// the force-constants matrix on
//
//		Eigen::Tensor<double,8> ukvec(18*numAtoms,
//				qCoarseGrid(0), qCoarseGrid(1), qCoarseGrid(2),
//				3, 3, numAtoms, numAtoms);
//		ukvec.setZero();
//
//		Eigen::Tensor<double,7> frc_new(qCoarseGrid(0), qCoarseGrid(1),
//				qCoarseGrid(2), 3, 3, numAtoms, numAtoms);
//		frc_new = forceConstants;
//
//		Eigen::array<long,8> offset8;
//		Eigen::array<long,8> extent8;
//
//		p = 0;
//		for ( int i=0; i<3; i++ ) {
//			for ( int j=0; j<3; j++ ) {
//				for ( int na=0; na<numAtoms; na++ ) {
//					// These are the 3*3*nat vectors associated with the
//					// translational acoustic sum rules
//
//					offset8 = {p,0,0,0,i,j,na,0};
//					extent8 = {p+1,qCoarseGrid(0), qCoarseGrid(1),
//							qCoarseGrid(2),i+1,j+1,na+1,numAtoms};
//					// ukvec(p,:,:,:,i,j,na,:) = 1.0;
//					ukvec.slice(offset8, extent8) = 1.;
//					p += 1;
//				}
//			}
//		}
//
//		//  if (n.eq.4) then
//		//     do i=1,3
//		//        do na=1,nat
//		//           ! These are the 3*nat vectors associated with the
//		//           ! single rotational sum rule (1D system)
//		//           p=p+1
//		//           do nb=1,nat
//		//              u(p) % vec (:,:,:,i,MOD(axis,3)+1,na,nb)=-tau(MOD(axis+1,3)+1,nb)
//		//              u(p) % vec (:,:,:,i,MOD(axis+1,3)+1,na,nb)=tau(MOD(axis,3)+1,nb)
//		//           enddo
//		//           !
//		//        enddo
//		//     enddo
//		//  endif
//
//		//  if (n.eq.6) then
//		//     do i=1,3
//		//        do j=1,3
//		//           do na=1,nat
//		//              ! These are the 3*3*nat vectors associated with the
//		//              ! three rotational sum rules (0D system - typ. molecule)
//		//              p=p+1
//		//              do nb=1,nat
//		//                 u(p) % vec (:,:,:,i,MOD(j,3)+1,na,nb)=-tau(MOD(j+1,3)+1,nb)
//		//                 u(p) % vec (:,:,:,i,MOD(j+1,3)+1,na,nb)=tau(MOD(j,3)+1,nb)
//		//              enddo
//		//              !
//		//           enddo
//		//        enddo
//		//     enddo
//		//  endif
//
//		Eigen::Tensor<int,3> ind_v(9*numAtoms*numAtoms*qCoarseGrid(0)*qCoarseGrid(1)
//				*qCoarseGrid(2),2,7);
//		Eigen::Tensor<double,2> v(9*numAtoms*numAtoms*qCoarseGrid(0)
//				*qCoarseGrid(1)*qCoarseGrid(2),2);
//
//		int m = 0;
//		int q, l;
//		for ( int n1=0; n1<qCoarseGrid(0); n1++ ) {
//			for ( int n2=0; n2<qCoarseGrid(1); n2++ ) {
//				for ( int n3=0; n3<qCoarseGrid(2); n3++ ) {
//					for ( int i=0; i<3; i++ ) {
//						for ( int j=0; j<3; j++ ) {
//							for ( int na=0; na<numAtoms; na++ ) {
//								for ( int nb=0; nb<numAtoms; nb++ ) {
//
//									// These are the vectors associated with
//									// the symmetry constraints
//									q = 1;
//									l = 0;
//
//									while ( ( l <= m ) && ( q != 0 ) ) {
//										if ( (ind_v(l,0,0)==n1) &&
//												(ind_v(l,0,1)==n2) &&
//												(ind_v(l,0,2)==n3) &&
//												(ind_v(l,0,3)==i) &&
//												(ind_v(l,0,4)==j) &&
//												(ind_v(l,0,5)==na) &&
//												(ind_v(l,0,6)==nb)) { q = 0; };
//										if ( (ind_v(l,1,0)==n1) &&
//												(ind_v(l,1,1)==n2) &&
//												(ind_v(l,1,2)==n3) &&
//												(ind_v(l,1,3)==i) &&
//												(ind_v(l,1,4)==j) &&
//												(ind_v(l,1,5)==na) &&
//												(ind_v(l,1,6)==nb)) { q = 0; };
//										l += 1;
//									}
//									if ( (n1==(qCoarseGrid(0)-n1) % qCoarseGrid(0)) &&
//											(n2==(qCoarseGrid(1)-n2) % qCoarseGrid(1)) &&
//											(n3==(qCoarseGrid(2)-n3) % qCoarseGrid(2)) &&
//											(i==j) && (na==nb)) { q = 0; };
//									if ( q != 0 ) {
//										ind_v(m,0,0) = n1;
//										ind_v(m,0,1) = n2;
//										ind_v(m,0,2) = n3;
//										ind_v(m,0,3) = i;
//										ind_v(m,0,4) = j;
//										ind_v(m,0,5) = na;
//										ind_v(m,0,6) = nb;
//										v(m,0) = 1.0 / sqrt(2.);
//										ind_v(m,1,0) = (qCoarseGrid(0)-n1) % qCoarseGrid(0);
//										ind_v(m,1,1) = (qCoarseGrid(1)-n2) % qCoarseGrid(1);
//										ind_v(m,1,2) = (qCoarseGrid(2)-n3) % qCoarseGrid(2);
//										ind_v(m,1,3) = j;
//										ind_v(m,1,4) = i;
//										ind_v(m,1,5) = nb;
//										ind_v(m,1,6) = na;
//										v(m,1) = - 1.0 / sqrt(2.);
//										m += 1; // for next iteration
//									}
//								}
//							}
//						}
//					}
//				}
//			}
//		}
//
//		// Gram-Schmidt orthonormalization of the set of vectors created.
//		// Note that the vectors corresponding to symmetry constraints are already
//		// orthonormalized by construction.
//
//		int n_less = 0;
//		Eigen::Tensor<double,7> w(qCoarseGrid(0),qCoarseGrid(1),qCoarseGrid(2),
//				3,3,numAtoms,numAtoms);
//		Eigen::Tensor<double,7> x(qCoarseGrid(0),qCoarseGrid(1),qCoarseGrid(2),
//				3,3,numAtoms,numAtoms);
//		Eigen::Tensor<double,7> x2(qCoarseGrid(0),qCoarseGrid(1),qCoarseGrid(2),
//				3,3,numAtoms,numAtoms);
//		Eigen::MatrixXi u_less(6*3*numAtoms);
//		u_less.setZero();
//
//		Eigen::Tensor<double,1> v_slice(2);
//		Eigen::Tensor<int,2> ind_v_slice(2,7);
//
//		Eigen::array<long,2> offset2;
//		Eigen::array<long,2> extent2;
//		Eigen::array<long,3> offset3;
//		Eigen::array<long,3> extent3;
//
//		int n1, n2, n3, i, j, na, nb, na1, i1, j1;
//
//		for ( int k=0; k<p; k++ ) {
//
//			offset8 = {k, 0, 0, 0, 0, 0, 0, 0};
//			extent8 = {k+1, qCoarseGrid(0), qCoarseGrid(1), qCoarseGrid(2),
//							3, 3, numAtoms, numAtoms};
//			w = ukvec.slice(offset8,extent8);
//			x = ukvec.slice(offset8,extent8);
//			// w = ukvec(k,:,:,:,:,:,:,:);
//			// x = ukvec(k,:,:,:,:,:,:,:);
//
//			for ( int l=0; l<m; l++ ) {
//
//				offset2 = {l,0};
//				extent2 = {l+1,2};
//				offset3 = {l,0,0};
//				extent3 = {l+1,2,7};
//				v_slice = v.slice(extent2,offset2);
//				ind_v_slice = ind_v.slice(extent3,offset3);
//
//				sp2(x, v_slice, ind_v_slice, scal);
//				for ( int r=0; r<2; r++ ) {
//					n1 = ind_v(l,r,0);
//					n2 = ind_v(l,r,1);
//					n3 = ind_v(l,r,2);
//					i = ind_v(l,r,3);
//					j = ind_v(l,r,4);
//					na = ind_v(l,r,5);
//					nb = ind_v(l,r,6);
//					w(n1,n2,n3,i,j,na,nb) -= scal * v(l,r);
//				}
//			}
//			if ( k+1 <= 9*numAtoms ) {
//				na1 = k % numAtoms;
//				if ( na1 == 0 ) na1 = numAtoms;
//				j1 = ((k-na1)/numAtoms) % 3;
//				i1 = ((((k-na1)/numAtoms)-j1)/3) % 3;
//			} else {
//				q = k - 9*numAtoms;
//				//        if (n.eq.4) then
//				//           na1=MOD(q,nat)
//				//           if (na1.eq.0) na1=nat
//				//           i1=MOD((q-na1)/nat,3)+1
//				//        else
//				na1 = q % numAtoms;
//				if ( na1 == 0 ) { na1 = numAtoms; }
//				j1 = ((q-na1)/numAtoms) % 3;
//				i1 = ((((q-na1)/numAtoms)-j1)/3) % 3;
//				//        endif
//			}
//			for ( int q=0; q<k-1; q++ ) {
//				r = 1;
//				for ( int i_less=0; i_less<n_less; i_less++ ) {
//					if ( u_less(i_less) == q ) {r = 0;};
//				}
//				if ( r != 0 ) {
//					offset8 = {q, 0, 0, 0, 0, 0, 0, 0};
//					extent8 = {q+1, qCoarseGrid(0), qCoarseGrid(1),
//							qCoarseGrid(2), 3, 3, numAtoms, numAtoms};
//					// ukvec(q,:,:,:,:,:,:,:)
//					x2 = ukvec.slice(offset8, extent8);
//					sp3(x, x2, i1, na1, scal);
//					w -= scal * x2;
//				}
//			}
//			sp1(w, w, norm2);
//			if ( norm2 > 1.0e-16 ) {
//				offset8 = {k, 0, 0, 0, 0, 0, 0, 0};
//				extent8 = {k+1, qCoarseGrid(0), qCoarseGrid(1), qCoarseGrid(2),
//								3, 3, numAtoms, numAtoms};
//				// ukvec(k,:,:,:,:,:,:,:) = w / sqrt(norm2);
//				ukvec.slice(offset8, extent8) = w / sqrt(norm2);
//			} else {
//				u_less(n_less) = k;
//				n_less += 1;
//			}
//		}
//
//		// Projection of the force-constants "vector" on the orthogonal of the
//		// subspace of the vectors verifying the sum rules and symmetry contraints
//
//
//		w.setZero();
//		for ( int l=0; l<m; l++ ) {
//			offset2 = {l,0};
//			extent2 = {l+1,2};
//			offset3 = {l,0,0};
//			extent3 = {l+1,2,7};
//			v_slice = v.slice(extent2,offset2);
//			ind_v_slice = ind_v.slice(extent3,offset3);
//			sp2(frc_new, v_slice, ind_v_slice, scal);
//			// sp2(frc_new, v(l,:), ind_v(l,:,:), scal);
//			for ( int r=0; r<2; r++ ) {
//				n1 = ind_v(l,r,0);
//				n2 = ind_v(l,r,1);
//				n3 = ind_v(l,r,2);
//				i = ind_v(l,r,3);
//				j = ind_v(l,r,4);
//				na = ind_v(l,r,5);
//				nb = ind_v(l,r,6);
//				w(n1,n2,n3,i,j,na,nb) += scal * v(l,r);
//			}
//		}
//
//		for ( int k=0; k<p; k++ ) {
//			r = 1;
//			for ( int i_less=0; i_less<n_less; i_less++ ) {
//				if ( u_less(i_less) == k ) { r = 0; }
//			}
//			if ( r != 0 ) {
//				// ukvec(k,:,:,:,:,:,:,:);
//				offset8 = {k, 0, 0, 0, 0, 0, 0, 0};
//				extent8 = {k+1, qCoarseGrid(0), qCoarseGrid(1), qCoarseGrid(2),
//								3, 3, numAtoms, numAtoms};
//				x = ukvec.slice(offset8,extent8);
//				sp1(x, frc_new, scal);
//				w += scal * ukvec.slice(offset8,extent8);
//			}
//		}
//
//		// Final substraction of the former projection to the initial frc, to get
//		// the new "projected" frc
//
//		frc_new -= w;
//		sp1(w, w, norm2);
//		std::cout << "Difference between old and new force-constants: "
//				<< sqrt(norm2) << std::endl;
//
//		forceConstants = frc_new;
//	}
//}


//void PhononH0::sp_zeu(Eigen::Tensor<double,3>& zeu_u,
//		Eigen::Tensor<double,3>& zeu_v,
//		double& scal) {
//	// does the scalar product of two effective charges matrices zeu_u and zeu_v
//	// (considered as vectors in the R^(3*3*nat) space, and coded in the usual way)
//
//	scal = 0.;
//	for ( int i=0; i<3; i++) {
//		for ( int j=0; j<3; j++) {
//			for ( int na=0; na<numAtoms; na++) {
//				scal += zeu_u(i,j,na) * zeu_v(i,j,na);
//			}
//		}
//	}
//}
//
//
//void PhononH0::sp1(Eigen::Tensor<double,7>& u,
//		Eigen::Tensor<double,7>& v,
//		double& scal) {
//	// does the scalar product of two force-constants matrices u and v
//	// (considered as vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space,
//	// and coded in the usual way)
//	scal = 0.;
//	for ( int n1=0; n1<qCoarseGrid(0); n1++ ) {
//		for ( int n2=0; n2<qCoarseGrid(1); n2++ ) {
//			for ( int n3=0; n3<qCoarseGrid(2); n3++ ) {
//				for ( int i=0; i<3; i++ ) {
//					for ( int j=0; j<3; j++ ) {
//						for ( int na=0; na<numAtoms; na++ ) {
//							for ( int nb=0; nb<numAtoms; nb++ ) {
//								scal += u(n1,n2,n3,i,j,na,nb)
//    	    									* v(n1,n2,n3,i,j,na,nb);
//							}
//						}
//					}
//				}
//			}
//		}
//	}
//}
//
//
//void PhononH0::sp2(Eigen::Tensor<double,7>u,
//		Eigen::Tensor<double,1>& v,
//		Eigen::Tensor<int,2>& ind_v, double& scal) {
//	// Does the scalar product of two force-constants matrices u and v
//	// (considered as vectors in the R^(3*3*nat*nat*nr1*nr2*nr3) space).
//	// u is coded in the usual way but v is coded as explained when defining
//	// the vectors corresponding to the symmetry constraints
//	scal = 0.;
//	for ( int i=0; i<2; i++ ) {
//		scal += u(ind_v(i,0),ind_v(i,1),ind_v(i,2),ind_v(i,3),ind_v(i,4),
//				ind_v(i,5), ind_v(i,6)) * v(i);
//	}
//}
//
//
//void PhononH0::sp3(Eigen::Tensor<double,7>& u,
//		Eigen::Tensor<double,7>& v,
//		const int& i, const int& na, double& scal) {
//	// like sp1, but in the particular case when u is one of the u(k)%vec
//	// defined in set_asr (before orthonormalization). In this case most of the
//	// terms are zero (the ones that are not are characterized by i and na), so
//	// that a lot of computer time can be saved (during Gram-Schmidt).
//	scal = 0.;
//	for ( int n1=0; n1<qCoarseGrid(0); n1++ ) {
//		for ( int n2=0; n2<qCoarseGrid(1); n2++ ) {
//			for ( int n3=0; n3<qCoarseGrid(2); n3++ ) {
//				for ( int j=0; j<3; j++ ) {
//					for ( int nb=0; nb<numAtoms; nb++ ) {
//						scal += u(n1,n2,n3,i,j,na,nb) * v(n1,n2,n3,i,j,na,nb);
//					}
//				}
//			}
//		}
//	}
//}





