#ifndef PHINTERACTION_H
#define PHINTERACTION_H

#include <complex>
#include "eigen.h"
#include "state.h"
#include "constants.h"
#include "crystal.h"
#include "utilities.h"

// * Class to calculate the coupling of a 3-phonon interaction given three
// * phonon modes.
class Interaction3Ph{
private:
	Crystal & crystal;
	long numTriplets;
	Eigen::Tensor<double,4> ifc3Tensor;
	Eigen::Tensor<double,3> cellPositions;
	Eigen::Tensor<long,2> displacedAtoms;

	template <typename A,typename B,typename C>
	std::tuple<Eigen::Tensor<double,3>, Eigen::Tensor<double,3>>
			calcCouplingSquared(A & state1, B & state2,
					C & state3Plus, C & state3Mins);

	// we set to zero rates for scattering with states of energies <0.001 cm^-1
	const double energyCutoff = 0.001 / ryToCmm1;
public:

	Interaction3Ph(Crystal & crystal_,
			long & numTriplets_,
			Eigen::Tensor<double,4> & ifc3Tensor_,
			Eigen::Tensor<double,3> & cellPositions_,
			Eigen::Tensor<long,2> & displacedAtoms_); // default constructor
	Interaction3Ph(const Interaction3Ph & that); // copy constructor
	Interaction3Ph & operator=(const Interaction3Ph & that);// assignment op

	// this is an interface: we can compute it on the fly or read the cache.
	template <typename A,typename B,typename C>
	std::tuple<Eigen::Tensor<double,3>, Eigen::Tensor<double,3>>
			getCouplingSquared(A & state1, B & state2,
					C & state3Plus, C & state3Mins);

	/**
	 * Calculate single three-phonon matrix element (V^{+}/V^{-1}).
	 *
	 * Method for calculating the matrix element of a single, plus
	 * or minus, three-phonon scattering process.
	 *
	 * @param[in] interactingPhonons phononTriplet data structure
	 *  containing the three phonon modes (in the full Brillouin zone)
	 *  taking part in an interaction event.
	 * @param[in] phononEigenvectors Eigenvectors of all phonons in the
	 *  full Brillouin zone.
	 * @param[in] numTriplets Number of triplets considered in the supercell
	 *  for the third order force constants calculation.
	 * @param[in] ifc3Tensor Third order force constants.
	 * @param[in] cellPositions Cartesian coordinates of the 2nd and 3rd unitcells.
	 * @param[in] displacedAtoms Index of the displaced atom for every triplet
	 *  of unitcells.
	 * @param[in] crystal Contains of crystal information
	 * @param[in] procType Character '+' or '-' to choose the +/- type process.
	 * @return The complex matrix element.
	 */

//	double calculateSingleV(const PhononTriplet &interactingPhonons, const Eigen::MatrixXd &q,
//			const int numTriplets, const Eigen::Tensor<double,4> &ifc3Tensor,
//			const Eigen::Tensor<double,3> &cellPositions,
//			const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo,
//			const char procType);

//  void calculateAllVminus(const int *grid, const PhononMode &mode,
//			    const Eigen::MatrixXd &qFBZ,
//			    const Eigen::Tensor<complex<double>,3> &ev, const int numTriplets,
//			    const Eigen::Tensor<double,4> &ifc3Tensor,
//			    const Eigen::Tensor<double,3> &cellPositions,
//			    const Eigen::Tensor<int,2> &displacedAtoms,const CrystalInfo &crysInfo);
//
//  void calculateAllW(const double T,const int *grid, const PhononMode &mode,
//		     const Eigen::MatrixXi &indexMesh, const CrystalInfo &crysInfo,
//		     const Eigen::MatrixXd omega, const TetraData tetra);
};


template <typename A,typename B,typename C>
std::tuple<Eigen::Tensor<double,3>, Eigen::Tensor<double,3>>
		Interaction3Ph::getCouplingSquared(A & state1, B & state2,
				C & state3Plus, C & state3Mins) {
	return calcCouplingSquared(state1, state2, state3Plus, state3Mins);
}

template <typename A,typename B,typename C>
std::tuple<Eigen::Tensor<double,3>, Eigen::Tensor<double,3>>
		Interaction3Ph::calcCouplingSquared(A & state1, B & state2,
				C & state3Plus, C & state3Mins ) {

	Eigen::Vector3d cell2Pos, cell3Pos;

	//Cartesian phonon wave vectors: q1,q2,q3
	auto q2 = state2.getCoords("cartesian");
	auto q3Plus = state3Plus.getCoords("cartesian");
	auto q3Mins = state3Mins.getCoords("cartesian");

	auto energies1 = state1.getEnergies();
	auto energies2 = state2.getEnergies();
	auto energies3Plus = state3Plus.getEnergies();
	auto energies3Mins = state3Mins.getEnergies();

	//phonon branches:
	// we allow the number of bands to be different in each direction
	long nb1 = energies1.size(); // <- numBands
	long nb2 = energies2.size();
	long nb3Plus = energies3Plus.size();
	long nb3Mins = energies3Mins.size();

	Eigen::Tensor<std::complex<double>,3> vPlus(nb1,nb2,nb3Plus);
	Eigen::Tensor<std::complex<double>,3> vMins(nb1,nb2,nb3Mins);
	vPlus.setZero();
	vMins.setZero();

	if ( false ) {

		// this is the first implementation of this function
		// it's reasonable, but it's about 30 times slower than the
		// implementation below

		Eigen::Tensor<std::complex<double>,3> ev1, ev2, ev3Plus, ev3Mins;

		state1.getEigenvectors(ev1); // size (3,numAtoms,numBands)
		state2.getEigenvectors(ev2);
		state3Plus.getEigenvectors(ev3Plus);
		state3Mins.getEigenvectors(ev3Mins);

		for ( long it=0; it<numTriplets; it++ ) { // sum over all triplets

			Eigen::Tensor<std::complex<double>,3> v0Plus(nb1,nb2,nb3Plus);
			Eigen::Tensor<std::complex<double>,3> v0Mins(nb1,nb2,nb3Mins);
			v0Plus.setZero();
			v0Mins.setZero();

			for ( int ib1=0; ib1<nb1; ib1++ ) {
				for ( int ib2=0; ib2<nb2; ib2++ ) {
					for ( int ic3 : {0,1,2} ) {
						for ( int ic2 : {0,1,2} ) {
							for ( int ic1 : {0,1,2} ) {
								// note: we should use ev2 without conjugation
								// but computed at -q2. Instead we use it at
								// +q2 and note that z^*(q)=z(-q)
								for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
									v0Plus(ib1,ib2,ib3) +=
											ifc3Tensor(it,ic1,ic2,ic3)
											* ev1(ic1,displacedAtoms(it,0),ib1)
									* std::conj(ev2(ic2,displacedAtoms(it,1),ib2))
									* std::conj(ev3Plus(ic3,displacedAtoms(it,2),ib3));
								}

								for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
									v0Mins(ib1,ib2,ib3) +=
											ifc3Tensor(it,ic1,ic2,ic3)
											* ev1(ic1,displacedAtoms(it,0),ib1)
									* std::conj(ev2(ic2,displacedAtoms(it,1),ib2))
									* std::conj(ev3Mins(ic3,displacedAtoms(it,2),ib3));
								}
							}
						}
					}
				}
			}

			for ( int ic : {0,1,2} ) {
				cell2Pos(ic) = cellPositions(it,0,ic);
				cell3Pos(ic) = cellPositions(it,1,ic);
			}

			// As a convention, the first primitive cell in the triplet is
			// restricted to the origin, so the phase for that cell is unity.

			double arg = - q2.dot(cell2Pos) - q3Plus.dot(cell3Pos);
			std::complex<double> phasePlus = exp( complexI * arg );
			arg = - q2.dot(cell2Pos) - q3Mins.dot(cell3Pos);
			std::complex<double> phaseMins = exp( complexI * arg );

			for ( int ib1=0; ib1<nb1; ib1++ ) {
				for ( int ib2=0; ib2<nb2; ib2++ ) {
					for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
						// case +
						vPlus(ib1,ib2,ib3) += v0Plus(ib1,ib2,ib3) * phasePlus;
					}
					for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
						// case -
						vMins(ib1,ib2,ib3) += v0Mins(ib1,ib2,ib3) * phaseMins;
					}
				}
			}
		}

	} else {
		// this is a bit more convoluted than the implementation above,
		// but it should be faster, as we break loops in two sections

		long numAtoms = crystal.getNumAtoms();
		long numBands = numAtoms * 3;

		Eigen::MatrixXcd ev1, ev2, ev3Plus, ev3Mins;

		state1.getEigenvectors(ev1); // size (3,numAtoms,numBands)
		state2.getEigenvectors(ev2);
		state3Plus.getEigenvectors(ev3Plus);
		state3Mins.getEigenvectors(ev3Mins);

		// note: tmp* is a tensor over cartesian and atomic indices
		// (whose size coincides with the band numbers)
		Eigen::Tensor<std::complex<double>,3> tmpPlus(numBands, numBands, numBands);
		Eigen::Tensor<std::complex<double>,3> tmpMins(numBands, numBands, numBands);
		tmpPlus.setZero();
		tmpMins.setZero();

		for ( long it=0; it<numTriplets; it++ ) { // sum over all triplets

			for ( int ic : {0,1,2} ) {
				cell2Pos(ic) = cellPositions(it,0,ic);
				cell3Pos(ic) = cellPositions(it,1,ic);
			}

			// As a convention, the first primitive cell in the triplet is
			// restricted to the origin, so the phase for that cell is unity.

			double arg = - q2.dot(cell2Pos) - q3Plus.dot(cell3Pos);
			std::complex<double> phasePlus = exp( complexI * arg );
			arg = - q2.dot(cell2Pos) - q3Mins.dot(cell3Pos);
			std::complex<double> phaseMins = exp( complexI * arg );

			long ind1, ind2, ind3;

			for ( int ic1 : {0,1,2} ) {
				ind1 = compress2Indeces(ic1, displacedAtoms(it,0), 3, numAtoms);
				for ( int ic2 : {0,1,2} ) {
					ind2 = compress2Indeces(ic2, displacedAtoms(it,1), 3, numAtoms);
					for ( int ic3 : {0,1,2} ) {
						ind3 = compress2Indeces(ic3, displacedAtoms(it,2), 3, numAtoms);
						tmpPlus(ind1, ind2, ind3) +=
								ifc3Tensor(it,ic1,ic2,ic3) * phasePlus;
						tmpMins(ind1, ind2, ind3) +=
								ifc3Tensor(it,ic1,ic2,ic3) * phaseMins;
					}
				}
			}
		}

		// this 6 loops are easier to read, but can be optimized as done below
		// and we break them into 3 loops with 4 nested loops
		// ( ~numBands^2 faster)
//		for ( int iac1=0; iac1<numBands; iac1++ ) {
//			for ( int iac2=0; iac2<numBands; iac2++ ) {
//				for ( int iac3=0; iac3<numBands; iac3++ ) {
//					for ( int ib1=0; ib1<nb1; ib1++ ) {
//						for ( int ib2=0; ib2<nb2; ib2++ ) {
//							for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
//								vPlus(ib1,ib2,ib3) +=
//										tmpPlus(iac1,iac2,iac3)
//										* ev1(iac1,ib1)
//										* std::conj(ev2(iac2,ib2))
//										* std::conj(ev3Plus(iac3,ib3));
//							}
//							for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
//								vMins(ib1,ib2,ib3) +=
//										tmpMins(iac1,iac2,iac3)
//										* ev1(iac1,ib1)
//										* std::conj(ev2(iac2,ib2))
//										* std::conj(ev3Mins(iac3,ib3));
//							}
//						}
//					}
//				}
//			}
//		}
		Eigen::Tensor<std::complex<double>,3> tmp1Plus(nb1,numBands,numBands);
		Eigen::Tensor<std::complex<double>,3> tmp1Mins(nb1,numBands,numBands);
		tmp1Plus.setZero();
		tmp1Mins.setZero();
		for ( int iac2=0; iac2<numBands; iac2++ ) {
			for ( int iac3=0; iac3<numBands; iac3++ ) {
				for ( int iac1=0; iac1<numBands; iac1++ ) {
					for ( int ib1=0; ib1<nb1; ib1++ ) {
						tmp1Plus(ib1,iac2,iac3) += tmpPlus(iac1,iac2,iac3)
								* ev1(iac1,ib1);
						tmp1Mins(ib1,iac2,iac3) += tmpMins(iac1,iac2,iac3)
								* ev1(iac1,ib1);
					}
				}
			}
		}
		Eigen::Tensor<std::complex<double>,3> tmp2Plus(nb1,nb2,numBands);
		Eigen::Tensor<std::complex<double>,3> tmp2Mins(nb1,nb2,numBands);
		tmp2Plus.setZero();
		tmp2Mins.setZero();
		for ( int ib1=0; ib1<nb1; ib1++ ) {
			for ( int iac3=0; iac3<numBands; iac3++ ) {
				for ( int iac2=0; iac2<numBands; iac2++ ) {
					for ( int ib2=0; ib2<nb2; ib2++ ) {
						tmp2Plus(ib1,ib2,iac3) += tmp1Plus(ib1,iac2,iac3)
								* std::conj(ev2(iac2,ib2));
						tmp2Mins(ib1,ib2,iac3) += tmp1Mins(ib1,iac2,iac3)
								* std::conj(ev2(iac2,ib2));
					}
				}
			}
		}
		for ( int ib1=0; ib1<nb1; ib1++ ) {
			for ( int ib2=0; ib2<nb2; ib2++ ) {
				for ( int iac3=0; iac3<numBands; iac3++ ) {
					for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
						vPlus(ib1,ib2,ib3) += tmp2Plus(ib1,ib2,iac3)
								* std::conj(ev3Plus(iac3,ib3));
					}
					for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
						vMins(ib1,ib2,ib3) += tmp2Mins(ib1,ib2,iac3)
								* std::conj(ev3Mins(iac3,ib3));
					}
				}
			}
		}
	}

	Eigen::Tensor<double,3> couplingPlus(nb1,nb2,nb3Plus);
	Eigen::Tensor<double,3> couplingMins(nb1,nb2,nb3Mins);
	couplingPlus.setZero();
	couplingMins.setZero();

	for ( int ib1=0; ib1<nb1; ib1++ ) {
		double omega1 = energies1(ib1);
		if ( omega1 < energyCutoff ) continue;

		for ( int ib2=0; ib2<nb2; ib2++ ) {
			double omega2 = energies2(ib2);
			if ( omega2 < energyCutoff ) continue;

			// case +
			for ( int ib3=0; ib3<nb3Plus; ib3++ ) {
				double omega3Plus = energies3Plus(ib3);
				if ( omega3Plus < energyCutoff ) continue;

				double freqsPlus = omega1 * omega2 * omega3Plus;
				couplingPlus(ib1,ib2,ib3) = std::norm(vPlus(ib1,ib2,ib3))
						/ freqsPlus;
			}

			// case -
			for ( int ib3=0; ib3<nb3Mins; ib3++ ) {
				double omega3Mins = energies3Mins(ib3);
				if ( omega3Mins < energyCutoff ) continue;

				double freqsMins = omega1 * omega2 * omega3Mins;
				couplingMins(ib1,ib2,ib3) = std::norm(vMins(ib1,ib2,ib3))
						/ freqsMins;
			}
		}
	}
	return {couplingPlus, couplingMins};
}

#endif
