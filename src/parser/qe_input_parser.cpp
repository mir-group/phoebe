#include "qe_input_parser.h"
#include "constants.h"
#include <math.h> // round()
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>  // to declare istringstream
#include <algorithm> // to use .remove_if
#include <stdlib.h> // abs()

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>
#include <Eigen/Eigenvalues>
#include <Eigen/Core>
#include <assert.h>

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw () {
		return "Error reading the file input parameter";
	}
};

void error(std::string errMessage, int errCode) {
	if ( errCode != 0 ) {
		std::cout << errMessage << "\n";
		assert(errCode != 0);
	}
}

void warning(std::string errMessage) {
	std::cout << errMessage;
}

double calcVolume(const Eigen::Matrix3d& directUnitCell, const double alat)
{
	Eigen::Vector3d a1 = directUnitCell.col(0);
	Eigen::Vector3d a2 = directUnitCell.col(1);
	Eigen::Vector3d a3 = directUnitCell.col(2);
	double volume;
	volume = abs( a1.dot(( a2.cross(a3) )) );
	volume+= abs( a2.dot(( a3.cross(a1) )) );
	volume+= abs( a3.dot(( a1.cross(a2) )) );
	volume *= alat * alat * alat / 3.;
	return volume;
}

Eigen::MatrixXd calcReciprocalCell(const Eigen::Matrix3d& directUnitCell)
{
	Eigen::Matrix3d reciprocalCell = directUnitCell.inverse().transpose();
	return reciprocalCell;
}

Eigen::MatrixXd wsinit(const Eigen::Matrix3d& unitCell) {
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

double wsweight(const Eigen::VectorXd& r, const Eigen::MatrixXd& rws) {
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

void rigidBulk(const int nr1, const int nr2, const int nr3, const int numAtoms,
		Eigen::Tensor<std::complex<double>,4>& dyn, const Eigen::VectorXd& q,
		const Eigen::MatrixXd& tau, const Eigen::MatrixXd& dielectricMatrix,
		const Eigen::Tensor<double,3>& bornCharges, const Eigen::MatrixXd& bg,
		const double volumeUnitCell, const double alat, const bool loto_2d,
		const int sign) {
	// compute the rigid-ion (long-range) term for q
	// The long-range term used here, to be added to or subtracted from the
	// dynamical matrices, is exactly the same of the formula introduced in:
	// X. Gonze et al, PRB 50. 13035 (1994) . Only the G-space term is
	// implemented: the Ewald parameter alpha must be large enough to
	// have negligible r-space contribution

	//   integer ::  nr1, nr2, nr3    !  FFT grid
	//   integer ::  numAtoms              ! number of atoms
	//   complex(DP) :: dyn(3,3,numAtoms,numAtoms) ! dynamical matrix
	//   real(DP) &
	//        q(3),           &! q-vector
	//        tau(3,numAtoms),     &! atomic positions
	//        dielectricMatrix(3,3),     &! dielectric constant tensor
	//        bornCharges(3,3,numAtoms),   &! effective charges tensor
	//        !at(3,3),        &! direct     lattice basis vectors
	//        bg(3,3),        &! reciprocal lattice basis vectors
	//        volumeUnitCell,          &! unit cell volume
	//        alat,           &! cell dimension units
	//        sign             ! sign=+/-1.0 ==> add/subtract rigid-ion term
	//   logical :: loto_2d ! 2D LOTO correction

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

	if ( nr1 == 1 ) {
		nr1x = 0;
	} else {
		nr1x = (int) ( sqrt(geg) /
				sqrt( bg.col(0).transpose() * bg.col(0) )) + 1;
	}
	if ( nr2 == 1) {
		nr2x = 0;
	} else {
		nr2x = (int) ( sqrt(geg) /
				sqrt( bg.col(1).transpose() * bg.col(1) )) + 1;
	}
	if (nr3 == 1) {
		nr3x = 0;
	} else {
		nr3x = (int) ( sqrt(geg) /
				sqrt( bg.col(2).transpose() * bg.col(2) )) + 1;
	}

	if ( abs(sign) != 1. ) {
		error("wrong value for sign", 1);
	}

	if ( loto_2d ) {
		fac = sign * e2 * 4. * pi / volumeUnitCell * 0.5 * alat / bg(3,3);
		reff.setZero();
		for ( int i=0; i<2; i++ ) {
			for ( int j=0; j<2; j++ ) {
				reff(i,j) = dielectricMatrix(i,j) * 0.5 * twoPi / bg(3,3); // (eps)*c/2 in 2pi/a units
			}
		}
		for ( int i=0; i<2; i++ ) {
			reff(i,i) = reff(i,i) - 0.5 * twoPi / bg(3,3); // (-1)*c/2 in 2pi/a units
		}
	} else {
		fac = sign * e2 * 4. * pi / volumeUnitCell;
	}

	Eigen::VectorXd g(3);
	std::complex<double> phase;
	for ( int m1=-nr1x; m1<=nr1x; m1++ ) {
		for ( int m2=-nr2x; m2<=nr2x; m2++ ) {
			for ( int m3=-nr3x; m3<=nr3x; m3++ ) {
				g(0) = m1*bg(0,0) + m2*bg(0,1) + m3*bg(0,2);
				g(1) = m1*bg(1,0) + m2*bg(1,1) + m3*bg(1,2);
				g(2) = m1*bg(2,0) + m2*bg(2,1) + m3*bg(2,2);

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
								arg = ( (tau.row(na)-tau.row(nb)) * g ).value();
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
							arg = ( (tau.row(na)-tau.row(nb)) * g ).value();
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

void nonAnal(const int numAtoms, const Eigen::VectorXi& atomicSpecies,
		const Eigen::MatrixXd& dielectricMatrix, const Eigen::VectorXd& q,
		const Eigen::Tensor<double,3>& bornCharges,
		const double& volumeUnitCell,
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
					dyn(i,j,it,jt) += 4. * pi * e2 * zag(i) * zbg(j) / qeq / volumeUnitCell;
				}
			}
		}
	}
}

void nonAnalIFC(const int numAtoms, const Eigen::VectorXi& atomicSpecies,
		const Eigen::MatrixXd& dielectricMatrix, const Eigen::VectorXd& q,
		const Eigen::Tensor<double,3>& bornCharges, const double volumeUnitCell,
		const int nr1, const int nr2, const int nr3,
		Eigen::Tensor<std::complex<double>, 4>& f_of_q)
{
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
					f_of_q(i,j,na,nb) = 4. * pi * e2 * zag(i) * zbg(j) / qeq
							/ volumeUnitCell / (nr1*nr2*nr3);
				}
			}
		}
	}
}

void frc_blk(Eigen::Tensor<std::complex<double>, 4>& dyn, const Eigen::VectorXd& q,
		const Eigen::MatrixXd& tau, int numAtoms, int nr1, int nr2, int nr3,
		const Eigen::Tensor<double, 7>& forceConstants, const Eigen::MatrixXd& at,
		const Eigen::MatrixXd& rws,
		Eigen::Tensor<std::complex<double>, 4>& f_of_q, const bool frozenPhonon)
{
	// calculates the dynamical matrix at q from the (short-range part of the)
	// force constants21

	Eigen::VectorXd r(3), r_ws(3);
	double total_weight, arg, weight;
	static bool first = true;

	int n1ForCache, n2ForCache, n3ForCache;

	const int nr1Big = 2 * nr1, nr2Big = 2 * nr2, nr3Big = 2 * nr3;

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
								r(i) = n1 * at(i,0) + n2 * at(i,1) + n3*at(i,2);
								r_ws(i) = r(i) + tau(na,i) - tau(nb,i);
								if ( frozenPhonon ) {
									r_ws(i) = r(i) + tau(nb,i) - tau(na,i);
								}
							}

							x = wsweight(r_ws, rws);

							wscache(n3ForCache, n2ForCache, n1ForCache, nb, na)
							= x;
							total_weight += x;
						}
					}
				}

				if ( abs( total_weight - nr1 * nr2 * nr3 ) > 1.0e-8 ) {
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

						// sum over r vectors in the supercell - very safe range!
						for ( int i=0; i<3; i++ ) {
							r(i) = n1 * at(i,0) + n2 * at(i,1) + n3 * at(i,2);
						}

						weight = wscache(n3ForCache,n2ForCache,n1ForCache,nb,na);
						if ( weight > 0. ) {

							// find vector corresponding to r in original cell

							m1 = ( n1 + 1 ) % nr1;
							if ( m1 <= 0 ) { m1 += nr1; };
							m2 = ( n2 + 1 ) % nr2;
							if ( m2 <= 0 ) { m2 += nr2; };
							m3 = ( n3 + 1 ) % nr3;
							if ( m3 <= 0 ) { m3 += nr3; };
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
											(forceConstants(m1,m2,m3,ipol,jpol,na,nb) +
													f_of_q(ipol,jpol,na,nb)) *
													phase * weight;
								}
							}

						}
					}
				}
			}
		}
	}
}

void setupMat(const Eigen::VectorXd& q,
		Eigen::Tensor<std::complex<double>, 4>& dyn, const int numAtoms,
		const Eigen::MatrixXd& bg,
		const Eigen::MatrixXd& tau,
		const Eigen::MatrixXd& at,
		const double volumeUnitCell,
		const bool loto_2d,
		const Eigen::MatrixXd& dielectricMatrix,
		const Eigen::Tensor<double,3>& bornCharges,
		const Eigen::Tensor<double,7>& forceConstants,
		const int nr1, const int nr2, const int nr3,
		const bool has_zstar,
		const Eigen::MatrixXd &rws,
		const bool na_ifc,
		Eigen::Tensor<std::complex<double>, 4>& f_of_q,
		const bool frozenPhonon,
		const double& alat)
{
	// compute the dynamical matrix (the analytic part only)

	// Set matrix to zero
	dyn.setZero();

	frc_blk(dyn, q, tau, numAtoms, nr1, nr2, nr3,
			forceConstants, at, rws, f_of_q, frozenPhonon);

	if ( has_zstar && !na_ifc ) {
		rigidBulk(nr1, nr2, nr3, numAtoms, dyn, q, tau,
				dielectricMatrix, bornCharges, bg, volumeUnitCell, alat, loto_2d, +1.);
	}

}

void dyndiag(const int numAtoms,
		const Eigen::VectorXd speciesMasses,
		const Eigen::VectorXi atomicSpecies,
		Eigen::Tensor<std::complex<double>,4>& dyn,
		Eigen::VectorXd& w2, Eigen::Tensor<std::complex<double>,3>& z)
{
	// diagonalise the dynamical matrix
	// On input:  speciesMasses = masses, in amu
	// On output: w2 = energies, z = displacements

	// fill the two-indices dynamical matrix

	int nat3 = 3 * numAtoms;
	int iType, jType;
	std::complex<double> cx;
	Eigen::MatrixXcd dyn2Tmp(nat3, nat3);
	Eigen::MatrixXcd dyn2(nat3, nat3);

	for (int iat = 0; iat<numAtoms; iat++) {
		for (int jat = 0; jat<numAtoms; jat++) {
			for (int ipol = 0; ipol<3; ipol++) {
				for (int jpol = 0; jpol<3; jpol++) {
					dyn2Tmp(iat*3 + ipol, jat*3 + jpol) = dyn(ipol,jpol,iat,jat);
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

	w2 = eigensolver.eigenvalues();


	Eigen::VectorXd freq(nat3);
	for ( int i=0; i<nat3; i++ ) {
		if ( w2(i) < 0 ) {
			freq(i) = sqrt(-w2(i));
		} else {
			freq(i) = sqrt(w2(i));
		}
	}
	std::cout << freq.transpose() << "\n";
	std::cout << freq.transpose() * ryToCmm1 << "\n";

	Eigen::MatrixXcd zTemp = eigensolver.eigenvectors();

	//  displacements are eigenvectors divided by sqrt(speciesMasses)

	for ( int iband=0; iband<nat3; iband++ ) {
		for ( int iat=0; iat<numAtoms; iat++ ) {
			iType = atomicSpecies(iat);
			for ( int ipol=0; ipol<3; ipol++ ) {
				z(ipol,iat,iband) = zTemp(iat*3 + ipol, iband) / sqrt(speciesMasses(iType));
			}
		}
	}
};


void cryst_to_cart(Eigen::VectorXd& vec, const Eigen::MatrixXd& trmat,
		const int iflag)
{
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

void latgen(const int ibrav, Eigen::VectorXd& celldm, Eigen::Matrix3d& unitCell)
{
	//  !     sets up the crystallographic vectors a1, a2, and a3.
	//  !
	//  !     ibrav is the structure index:
	//  !       1  cubic P (sc)                8  orthorhombic P
	//  !       2  cubic F (fcc)               9  1-face (C) centered orthorhombic
	//  !       3  cubic I (bcc)              10  all face centered orthorhombic
	//  !       4  hexagonal and trigonal P   11  body centered orthorhombic
	//  !       5  trigonal R, 3-fold axis c  12  monoclinic P (unique axis: c)
	//  !       6  tetragonal P (st)          13  one face (base) centered monoclinic
	//  !       7  tetragonal I (bct)         14  triclinic P
	//  !     Also accepted:
	//  !       0  "free" structure          -12  monoclinic P (unique axis: b)
	//  !      -3  cubic bcc with a more symmetric choice of axis
	//  !      -5  trigonal R, threefold axis along (111)
	//  !      -9  alternate description for base centered orthorhombic
	//  !     -13  one face (base) centered monoclinic (unique axis: b)
	//  !      91  1-face (A) centered orthorombic
	//  !
	//  !     celldm are parameters which fix the shape of the unit cell
	//  !     volumeUnitCell is the unit-cell volume
	//  !
	//  !     NOTA BENE: all axis sets are right-handed
	//  !     Boxes for US PPs do not work properly with left-handed axis

	const double sr2 = 1.414213562373, sr3 = 1.732050807569;

	//  user-supplied lattice vectors

	Eigen::Vector3d a1, a2, a3;

	a1 = unitCell.col(0);
	a2 = unitCell.col(1);
	a3 = unitCell.col(2);

	if ( ibrav == 0 ) {
		if ( sqrt( a1.transpose()*a1 ) == 0 || sqrt( a2.transpose()*a2 ) == 0
				|| sqrt( a3.transpose()*a3 ) == 0 ) {
			error("wrong at for ibrav=0", 1);
		}
		if ( celldm(0) != 0. ) {
			// ... input at are in units of alat => convert them to a.u.
			unitCell *= celldm(0);
		} else {
			// ... input at are in atomic units: define celldm(1) from a1
			celldm(0) = sqrt( a1.transpose() * a1 );
		}
	} else {
		a1.setZero();
		a2.setZero();
		a3.setZero();
	}

	if ( celldm(0) <= 0. ) {
		error("wrong celldm(1)", 1 );
	}

	//  index of bravais lattice supplied

	if ( ibrav == 1 ) { // simple cubic lattice
		a1(0) = celldm(0);
		a2(1) = celldm(0);
		a3(2) = celldm(0);
	} else if (ibrav == 2) { //     fcc lattice
		double term = celldm(0) / 2.;
		a1(0) =-term;
		a1(2) = term;
		a2(1) = term;
		a2(2) = term;
		a3(0) =-term;
		a3(1) = term;
	} else if (abs(ibrav) == 3) { // bcc lattice
		double term = celldm(0) / 2.;
		for ( int ir=0; ir<3; ir++ ) {
			a1(ir) = term;
			a2(ir) = term;
			a3(ir) = term;
		} if ( ibrav < 0 ) {
			a1(0) = -a1(0);
			a2(1) = -a2(1);
			a3(2) = -a3(2);
		} else {
			a2(0) = -a2(0);
			a3(0) = -a3(0);
			a3(1) = -a3(1);
		}
	} else if ( ibrav == 4 ) {// hexagonal lattice
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		double cbya  = celldm(2);
		a1(1) = celldm(0);
		a2(1) =-celldm(0) / 2.;
		a2(2) = celldm(0) * sr3 / 2.;
		a3(3) = celldm(0) * cbya;

	} else if (abs(ibrav) == 5) { // trigonal lattice
		if ( celldm(3) <= -0.5 || celldm(3) >= 1. ) {
			error("wrong celldm(4)", abs(ibrav));
		}

		double term1 = sqrt(1. + 2. * celldm(3) );
		double term2 = sqrt(1. - celldm(3) );

		if ( ibrav == 5 ) { // threefold axis along c (001)
			a2(1) = sr2 * celldm(0) * term2 / sr3;
			a2(2) = celldm(0) * term1 / sr3;
			a1(0) = celldm(0) * term2 / sr2;
			a1(1) =-a1(0) / sr3;
			a1(2) = a2(2);
			a3(0) =-a1(0);
			a3(1) = a1(1);
			a3(2) = a2(2);
		} else if ( ibrav == -5 ) { // threefold axis along (111)
			// Notice that in the cubic limit (alpha=90, celldm(4)=0, term1=term2=1)
			// does not yield the x,y,z axis, but an equivalent rotated triplet:
			//   a/3 (-1,2,2), a/3 (2,-1,2), a/3 (2,2,-1)
			// If you prefer the x,y,z axis as cubic limit, you should modify the
			// definitions of a1(1) and a1(2) as follows:'
			// a1(1) = celldm(1)*(term1+2.0_dp*term2)/3.0_dp
			// a1(2) = celldm(1)*(term1-term2)/3.0_dp
			// (info by G. Pizzi and A. Cepellotti)
			a1(0) = celldm(0) * ( term1 - 2. * term2 ) / 3.;
			a1(1) = celldm(0) * ( term1 + term2 ) / 3.;
			a1(2) = a1(1);
			a2(0) = a1(2);
			a2(1) = a1(0);
			a2(2) = a1(1);
			a3(0) = a1(1);
			a3(1) = a1(2);
			a3(2) = a1(0);
		}
	} else if (ibrav == 6) { // tetragonal lattice
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", 6);
		}
		double cbya = celldm(2);
		a1(0) = celldm(0);
		a2(1) = celldm(0);
		a3(2) = celldm(0) * cbya;

	} else if (ibrav == 7) { // body centered tetragonal lattice
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", 7);
		}
		double cbya = celldm(2);
		a2(0) = celldm(0) / 2.;
		a2(1) = a2(0);
		a2(2) = cbya * celldm(0) / 2.;
		a1(0) = a2(0);
		a1(1) = - a2(0);
		a1(2) = a2(2);
		a3(0) = - a2(0);
		a3(1) = - a2(0);
		a3(2) = a2(2);
	} else if ( ibrav == 8 ) { // Simple orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. )
		{
			error("wrong celldm(3)", ibrav);
		}
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1);
		a3(2) = celldm(0) * celldm(2);
	} else if ( abs(ibrav) == 9) { // One face (base) centered orthorhombic lattice  (C type)
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", abs(ibrav));
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", abs(ibrav));
		}
		if ( ibrav == 9 ) {// old PWscf description
			a1(0) = 0.5 * celldm(0);
			a1(1) = a1(0) * celldm(1);
			a2(0) = - a1(0);
			a2(1) = a1(1);
		} else {// alternate description
			a1(0) =  0.5 * celldm(0);
			a1(1) = -a1(0) * celldm(1);
			a2(0) =  a1(0);
			a2(1) = -a1(1);
		}
		a3(2) = celldm(0) * celldm(2);
	} else if ( ibrav == 91 ) { // One face(base)centered orthorhombic lattice (A type)
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1) * 0.5;
		a2(2) = - celldm(0) * celldm(2) * 0.5;
		a3(1) = a2(1);
		a3(2) = - a2(2);
	} else if (ibrav == 10) {// All face centered orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		a2(0) = 0.5 * celldm(0);
		a2(1) = a2(0) * celldm(1);
		a1(0) = a2(0);
		a1(2) = a2(0) * celldm(2);
		a3(1) = a2(0) * celldm(1);
		a3(2) = a1(2);
	} else if (ibrav == 11) { // Body centered orthorhombic lattice
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		a1(0) = 0.5 * celldm(0);
		a1(1) = a1(0) * celldm(1);
		a1(2) = a1(0) * celldm(2);
		a2(0) = - a1(0);
		a2(1) = a1(1);
		a2(2) = a1(2);
		a3(0) = - a1(0);
		a3(1) = - a1(1);
		a3(2) = a1(2);
	} else if (ibrav == 12) { // Simple monoclinic lattice, unique (i.e. orthogonal to a) axis: c
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			error("wrong celldm(4)", ibrav);
		}
		double sen = sqrt( 1. - celldm(3)*celldm(3) );
		a1(0) = celldm(0);
		a2(0) = celldm(0) * celldm(1) * celldm(3);
		a2(1) = celldm(0) * celldm(1) * sen;
		a3(2) = celldm(0) * celldm(2);
	} else if ( ibrav == - 12 ) { // Simple monoclinic lattice, unique axis: b (more common)
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)",-ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)",-ibrav);
		}
		if ( abs(celldm(4))>=1. ) {
			error("wrong celldm(5)",-ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = celldm(0);
		a2(1) = celldm(0) * celldm(1);
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(2) = celldm(0) * celldm(2) * sen;
	} else if ( ibrav == 13 ) { // One face centered monoclinic lattice unique axis c
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			error("wrong celldm(4)", ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = 0.5 * celldm(0);
		a1(2) =-a1(0) * celldm(2);
		a2(0) = celldm(0) * celldm(1) * celldm(2);
		a2(1) = celldm(0) * celldm(1) * sen;
		a3(0) = a1(0);
		a3(2) =-a1(2);
	} else if ( ibrav == -13 ) { // One face centered monoclinic lattice unique axis b
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", -ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", -ibrav);
		}
		if ( abs(celldm(4)) >= 1. ) {
			error("wrong celldm(5)", -ibrav);
		}
		double sen = sqrt( 1. - celldm(4)*celldm(4) );
		a1(0) = 0.5 * celldm(0);
		a1(1) =-a1(0) * celldm(1);
		a2(0) = a1(0);
		a2(1) =-a1(1);
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(2) = celldm(0) * celldm(2) * sen;
	} else if (ibrav == 14) { // Triclinic lattice
		if ( celldm(1) <= 0. ) {
			error("wrong celldm(2)", ibrav);
		}
		if ( celldm(2) <= 0. ) {
			error("wrong celldm(3)", ibrav);
		}
		if ( abs(celldm(3)) >= 1. ) {
			error("wrong celldm(4)", ibrav);
		}
		if ( abs(celldm(4)) >= 1. ) {
			error("wrong celldm(5)", ibrav);
		}
		if ( abs(celldm(5)) >= 1. ) {
			error("wrong celldm(6)", ibrav);
		}
		double singam = sqrt( 1. - celldm(5)*celldm(5) );
		double term = ( 1. + 2. * celldm(3)*celldm(4)*celldm(5)
				- celldm(3)*celldm(3) - celldm(4)*celldm(4) - celldm(5)*celldm(5));
		if ( term < 0. )
		{
			error("celldm does not make sense, check your data", ibrav);
		}
		term = sqrt( term / ( 1. - celldm(5)*celldm(5) ) );
		a1(0) = celldm(0);
		a2(0) = celldm(0) * celldm(1) * celldm(5);
		a2(1) = celldm(0) * celldm(1) * singam;
		a3(0) = celldm(0) * celldm(2) * celldm(4);
		a3(1) = celldm(0) * celldm(2) * (celldm(3)-celldm(4)*celldm(5))/singam;
		a3(2) = celldm(0) * celldm(2) * term;

	} else {
		error("nonexistent bravais lattice", ibrav);
	}

	if ( ibrav != 0 ) {
		unitCell.col(0) = a1;
		unitCell.col(1) = a2;
		unitCell.col(2) = a3;
	}
}

void diagonalize(const Eigen::VectorXd& q, const int numAtoms,
		const Eigen::MatrixXd& dielectricMatrix,
		const Eigen::Tensor<double, 3>& bornCharges, const bool na_ifc,
		const Eigen::MatrixXi& atomicSpecies, // atomic types for each atom of the original cell
		const double volumeUnitCell,
		const int nr1, const int nr2, const int nr3,
		const Eigen::MatrixXd& at, const Eigen::MatrixXd& bg,
		const bool loto_2d,
		const Eigen::MatrixXd& tau,
		const Eigen::Tensor<double, 7> forceConstants,
		const bool has_zstar,
		const double& alat,
		const Eigen::MatrixXd& rws,
		const bool frozenPhonon,
		const Eigen::VectorXd speciesMasses
)
// to be executed at every q-point, to get phonon frequencies and wavevectors
{
	Eigen::Tensor<std::complex<double>, 4> dyn(3,3,numAtoms,numAtoms);
	dyn.setZero();
	//	Eigen::Tensor<std::complex<double>, 4> dyn_blk(3,3,numAtoms,numAtoms);
	//	dyn_blk.setZero();

	Eigen::Tensor<std::complex<double>,4> f_of_q(3,3,numAtoms,numAtoms);
	f_of_q.setZero();

	Eigen::VectorXd qhat(3);
	double qq;

	if ( na_ifc ) {
		qq = sqrt( q.transpose()*q ); // q is the qpoint coordinate
		if ( abs(qq) < 1.0e-8 ) {
			qq = 1.;
		}

		qhat = q / qq;

		nonAnalIFC(numAtoms, atomicSpecies, dielectricMatrix, qhat, bornCharges, volumeUnitCell,
				nr1, nr2, nr3, f_of_q);
	}

	setupMat(q, dyn, numAtoms, bg, tau,
			at, volumeUnitCell, loto_2d, dielectricMatrix, bornCharges,
			forceConstants, nr1,nr2,nr3, has_zstar, rws, na_ifc, f_of_q,
			frozenPhonon, alat);

	if ( !loto_2d && na_ifc ) {
		qhat = q.transpose() * at;
		if ( abs( qhat(0) - round(qhat(0) ) ) <= 1.0e-6 &&
				abs( qhat(1) - round(qhat(1) ) ) <= 1.0e-6 &&
				abs( qhat(2) - round(qhat(2) ) ) <= 1.0e-6 ) {
			// q = 0 : we need the direction q => 0 for the non-analytic part

			qq = sqrt( ( qhat.transpose()*qhat ).value() );
			if (qq != 0. ) {
				qhat /= qq;
			}
			nonAnal(numAtoms, atomicSpecies, dielectricMatrix, qhat, bornCharges, volumeUnitCell, dyn);
		}
	}

	Eigen::VectorXd w2(numAtoms*3);
	Eigen::Tensor<std::complex<double>,3> z(3,numAtoms,numAtoms*3);

	dyndiag(numAtoms, speciesMasses, atomicSpecies, dyn, w2, z);
}

std::vector<std::string> split(const std::string& s, char delimiter)
										{
	std::vector<std::string> tokens;
	std::string token;
	std::istringstream tokenStream(s);

	if ( delimiter == ' ' ) {
		for (std::string s; tokenStream >> s; ) {
			tokens.push_back(s);
		}
	} else {
		while (std::getline(tokenStream, token, delimiter)) {
			token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
					token.end());
			tokens.push_back(token);
		}
	}

	return tokens;
										}

void QEParser::parsePhHarmonic(std::string fileName) {
	//  Here we read the dynamical matrix of interatomic force constants
	//	in real space.
	//	Since the file is typically small, we don't worry about memory management

	std::string line;
	std::vector<std::string> lineSplit;

	// open input file
	std::ifstream infile(fileName);

	//  First line contains ibrav, celldm and other variables

	std::getline(infile, line);
	lineSplit = split(line, ' ');

	int numElements = std::stoi(lineSplit[0]);
	int numAtoms = std::stoi(lineSplit[1]);
	int ibrav = std::stoi(lineSplit[2]);

	Eigen::VectorXd celldm(6);
	celldm(0) = std::stod(lineSplit[3]);
	celldm(1) = std::stod(lineSplit[4]);
	celldm(2) = std::stod(lineSplit[5]);
	celldm(3) = std::stod(lineSplit[6]);
	celldm(4) = std::stod(lineSplit[7]);
	celldm(5) = std::stod(lineSplit[8]);

	Eigen::Matrix3d directUnitCell(3,3);
	if ( ibrav == 0 ) {
		//    	In this case, unitCell is written in the file, in angstroms
		for ( int i=0; i<3; i++ ) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				directUnitCell(i,j) = std::stod(lineSplit[j]); // / distanceRyToAng;
			}
		};
	};

	// generate the unit cell vectors (also for ibrav != 0)
	latgen(ibrav, celldm, directUnitCell);

	double alat = celldm(0);
	directUnitCell /= alat; // bring unit cell in units of the lattice parameter

	//  Next, we read the atomic species
	std::vector<std::string> speciesNames;
	Eigen::VectorXd speciesMasses(numElements);
	for ( int i=0; i<numElements; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, '\'');
		speciesNames.push_back(lineSplit[1]);
		speciesMasses(i) = std::stod(lineSplit[2]); // in rydbergs
	};

	//  we read the atomic positions
	Eigen::MatrixXd atomicPositions(numAtoms,3);
	Eigen::VectorXi atomicSpecies(numAtoms);
	for ( int i=0; i<numAtoms; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, ' ');
		atomicSpecies(i) = std::stoi(lineSplit[1]) - 1;
		atomicPositions(i,0) = std::stod(lineSplit[2]);
		atomicPositions(i,1) = std::stod(lineSplit[3]);
		atomicPositions(i,2) = std::stod(lineSplit[4]);
	}

	//  Read if hasDielectric
	std::getline(infile, line);
	line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
	bool hasDielectric;
	if ( line == "T" ) {
		hasDielectric = true;
	} else {
		hasDielectric = false;
	}

	//	if there are the dielectric info, we can read dielectric matrix
	//	and the Born charges
	Eigen::MatrixXd dielectricMatrix(3,3);
	dielectricMatrix.setZero();
	Eigen::Tensor<double,3> bornCharges(numAtoms, 3, 3);
	bornCharges.setZero();

	if ( hasDielectric ) {
		for ( int i=0; i<3; i++) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				dielectricMatrix(i,j) = std::stod(lineSplit[j]);
			}
		}

		for ( int iAtom=0; iAtom < numAtoms; iAtom++ ) {
			std::getline(infile, line);
			for ( int i=0; i<3; i++ ) {
				std::getline(infile, line);
				lineSplit = split(line, ' ');
				for ( int j=0; j<3; j++ ) {
					bornCharges(iAtom,i,j) = std::stod(lineSplit[j]);
				}
			}
		}
	}

	//	Now we parse the coarse q grid
	std::getline(infile, line);
	lineSplit = split(line, ' ');
	Eigen::VectorXi qCoarseGrid(3);
	qCoarseGrid(0) = std::stoi(lineSplit[0]);
	qCoarseGrid(1) = std::stoi(lineSplit[1]);
	qCoarseGrid(2) = std::stoi(lineSplit[2]);

	std::cout << qCoarseGrid[0];
	std::cout << qCoarseGrid[1];
	std::cout << qCoarseGrid[2];
	std::cout << numAtoms << "\n";

	Eigen::Tensor<double, 7> forceConstants(qCoarseGrid[0], qCoarseGrid[1],
			qCoarseGrid[2], 3, 3, numAtoms, numAtoms);

	int m1Test;
	int m2Test;
	int m3Test;
	double x;

	for ( int ic=0; ic<3; ic++ ) {
		for ( int jc=0; jc<3; jc++ ) {
			for ( int iat=0; iat<numAtoms; iat++ ) {
				for ( int jat=0; jat<numAtoms; jat++ ) {
					// a line containing ic, jc, iat, jat
					std::getline(infile, line);

					for ( int r3=0; r3<qCoarseGrid[2]; r3++ ) {
						for ( int r2=0; r2<qCoarseGrid[1]; r2++ ) {
							for ( int r1=0; r1<qCoarseGrid[0]; r1++ ) {
								std::getline(infile, line);
								istringstream iss(line);
								iss >> m1Test >> m2Test >> m3Test >> x;
								forceConstants(r1, r2, r3, ic, jc, iat, jat) = x;
							}
						}
					}
				}
			}
		}
	}

	infile.close();

	// Now we do postprocessing

	double volumeUnitCell = calcVolume(directUnitCell, alat);
	Eigen::Matrix3d reciprocalUnitCell = calcReciprocalCell(directUnitCell);

	if ( qCoarseGrid(0) <= 0 || qCoarseGrid(1) <= 0 || qCoarseGrid(2) <= 0 ) {
		error("qCoarseGrid smaller than zero", 1);
	}


	//	Now, let's try to diagonalize some points, and start debugging at q=0

	Eigen::VectorXd q(3);
	q << 0., 0., 0.2;

	bool na_ifc = false;
	bool loto_2d = false;
	bool frozenPhonon = false;

	//	bool has_zstar=false;  <- hasDielectric

	Eigen::Matrix3d directUnitCellSup(3,3);
	directUnitCellSup.col(0) = directUnitCell.col(0) * qCoarseGrid(0);
	directUnitCellSup.col(1) = directUnitCell.col(1) * qCoarseGrid(1);
	directUnitCellSup.col(2) = directUnitCell.col(2) * qCoarseGrid(2);

	Eigen::MatrixXd rws = wsinit(directUnitCellSup);

	diagonalize(q, numAtoms, dielectricMatrix,
			bornCharges, na_ifc,
			atomicSpecies,
			volumeUnitCell,
			qCoarseGrid(0), qCoarseGrid(1), qCoarseGrid(2),
			directUnitCell, reciprocalUnitCell,
			loto_2d,
			atomicPositions,
			forceConstants,
			hasDielectric,
			alat,
			rws,
			frozenPhonon,
			speciesMasses);

	return;
};


