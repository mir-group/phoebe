#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <math.h>
#include <iostream>
#include "phononH0.h"
#include "constants.h"

void error(std::string errMessage, int errCode) {
	if ( errCode != 0 ) {
		std::cout << errMessage << std::endl;
		assert(errCode != 0);
	}
}

void warning(std::string errMessage) {
	std::cout << errMessage;
}

Eigen::MatrixXd PhononH0::wsinit(const Eigen::Matrix3d& unitCell) {
	const int nx=2;
	int index = 0;
	const int nrwsx = 200;

	Eigen::MatrixXd tmpResult(3,nrwsx);

	for ( int ir=-nx; ir<=nx; ir++ ) {
		for ( int jr=-nx; jr<=nx; jr++ ) {
			for ( int kr=-nx; kr<=nx; kr++ ) {
				for ( int i=0; i<3; i++ ) {
					tmpResult(i,index) = unitCell(i,0)*ir + unitCell(i,1)*jr + unitCell(i,2)*kr;
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







