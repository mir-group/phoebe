#include "ph_scattering.h"
#include "constants.h"
#include "io.h"

PhScatteringMatrix::PhScatteringMatrix(Context & context_,
			StatisticsSweep & statisticsSweep_,
			FullBandStructure<FullPoints> & innerBandStructure_,
			FullBandStructure<FullPoints> & outerBandStructure_,
			Interaction3Ph * coupling3Ph_,
			PhononH0 * h0_) :
			ScatteringMatrix(context_, statisticsSweep_,
					innerBandStructure_, outerBandStructure_),
			coupling3Ph(coupling3Ph_), h0(h0_) {
//	couplingIsotope = couplingIsotope_;
//	couplingBoundary = couplingBoundary_;
	if ( &innerBandStructure != &outerBandStructure && h0 == nullptr ) {
		Error e("PhScatteringMatrix needs h0 for incommensurate grids");
	}
}

PhScatteringMatrix::PhScatteringMatrix(const PhScatteringMatrix & that) :
		ScatteringMatrix(that), coupling3Ph(that.coupling3Ph), h0(that.h0) {
//		couplingIsotope(that.couplingIsotope);
//		couplingBoundary(that.couplingBoundary);
}

PhScatteringMatrix & PhScatteringMatrix::operator=(
		const PhScatteringMatrix & that) {
	ScatteringMatrix::operator=(that);
	if ( this != &that ) {
		coupling3Ph = that.coupling3Ph;
//		couplingIsotope = that.couplingIsotope;
//		couplingBoundary = that.couplingBoundary;
		h0 = that.h0;
	}
	return *this;
}

// 3 cases:
// theMatrix and linedith is passed: we compute and store in memory the scatt
//       matrix and the diagonal
// inPopulation+outPopulation is passed: we compute the action of the
//       scattering matrix on the in vector, returning outVec = sMatrix*vector
// only linewidth is passed: we compute only the linewidths
void PhScatteringMatrix::builder(
		Eigen::MatrixXd & matrix, VectorBTE * linewidth,
		VectorBTE * inPopulation, VectorBTE * outPopulation) {

	// notes: + process is (1+2) -> 3
	//        - processes are (1+3)->2 and (3+2)->1

	const double energyCutoff = 0.001 / ryToCmm1; // discard states with small
	// energies (smaller than 0.001 cm^-1

	int switchCase = 0;
	if ( matrix.rows() != 0 && linewidth != nullptr &&
		inPopulation == nullptr && outPopulation == nullptr  ) {
		switchCase = 0;
	} else if ( matrix.rows() == 0 && linewidth == nullptr &&
		inPopulation != nullptr && outPopulation != nullptr  ) {
		switchCase = 1;
	} else if ( matrix.rows() == 0 && linewidth != nullptr &&
		inPopulation == nullptr && outPopulation == nullptr  ) {
		switchCase = 2;
	} else {
		Error e("builder3Ph found a non-supported case");
	}

	if ( (linewidth != nullptr) && (linewidth->dimensionality!= 1) ) {
		Error e("The linewidths shouldn't have dimensionality");
	}

	bool dontComputeQ3 = &innerBandStructure == &outerBandStructure;

	auto statistics = outerBandStructure.getStatistics();

	long numAtoms = innerBandStructure.getPoints().getCrystal().getNumAtoms();
	long numCalcs = statisticsSweep.getNumCalcs();

	long outerNumPoints = outerBandStructure.getNumPoints();
	long innerNumPoints = innerBandStructure.getNumPoints();

	// precompute Bose populations
	VectorBTE outerBose(context, outerBandStructure, 1);
	for ( auto ik=0; ik<outerNumPoints; ik++ ) {
		auto state = outerBandStructure.getState(ik);
		auto energies = state.getEnergies();
		for ( auto ib=0; ib<energies.size(); ib++ ) {
			long is = outerBandStructure.getIndex(ik,ib);
			auto energy = energies(ib);
			for ( long iCalc=0; iCalc<statisticsSweep.getNumCalcs(); iCalc++){
				double temperature = statisticsSweep.getCalcStatistics(iCalc
						).temperature;
				outerBose.data(iCalc,is) = statistics.getPopulation(energy,
						temperature);
			}
		}
	}
	VectorBTE innerBose(context, outerBandStructure, 1);
	if ( &innerBandStructure == &outerBandStructure ) {
		innerBose = outerBose;
	} else {
		for ( auto ik=0; ik<innerNumPoints; ik++ ) {
			auto state = innerBandStructure.getState(ik);
			auto energies = state.getEnergies();
			for ( auto ib=0; ib<energies.size(); ib++ ) {
				long is = innerBandStructure.getIndex(ik,ib);
				auto energy = energies(ib);
				for ( long iCalc=0; iCalc<statisticsSweep.getNumCalcs();
						iCalc++ ) {
					double temperature = statisticsSweep.getCalcStatistics(
							iCalc).temperature;
					outerBose.data(iCalc,is) = statistics.getPopulation(energy,
							temperature);
				}
			}
		}
	}

	// note: these variables are only needed in the loop
	// but since it's an expensive loop, we define them here once and for all
	long nb3Plus, nb3Mins, ind1, ind2;
	double ratePlus, rateMins1, rateMins2, deltaPlus, deltaMins1, deltaMins2;
	double en1, en2, en3Plus, en3Mins, bose1, bose2, bose3Plus, bose3Mins;
	Eigen::Tensor<double,3> couplingPlus, couplingMins;
	Eigen::VectorXd state3PlusEnergies, state3MinsEnergies;
	Eigen::Vector3d v1, v2;
	Eigen::MatrixXd bose3PlusData, bose3MinsData;
	Eigen::VectorXd eigvals3Plus, eigvals3Mins;
	Eigen::MatrixXcd eigvecs3Plus, eigvecs3Mins;
//	DetachedState states3Plus, states3Mins;

//	double * eigvals3Plus_, eigvals3Mins_;
//	std::complex<double> * eigvecs3Plus_, eigvecs3Mins_;
//	std::complex<double> * velPtr=nullptr;
//	std::vector<double> q3Plusv(3), q3Minsv(3);

	LoopPrint loopPrint("computing scattering matrix", "q-points",
			outerNumPoints);

	for( long iq1=0; iq1<outerNumPoints; iq1++ ) {
		loopPrint.update();

		// note: for computing linewidths on a path, we must distinguish
		// that q1 and q2 are on different meshes, and that q3+/- may not fall
		// into known meshes and therefore needs to be computed

		auto states1 = outerBandStructure.getState(iq1);
		auto q1 = states1.getPoint();
		auto state1Energies = states1.getEnergies();
		auto nb1 = state1Energies.size();

		for( long iq2=0; iq2<innerNumPoints; iq2++ ) {
			auto q2 = innerBandStructure.getPoint(iq2);
//			auto states2 = innerBandStructure.getState(q2);

			auto iq2Inv = innerBandStructure.getPoints().getIndexInverted(iq2);
			auto q2Inv = innerBandStructure.getPoint(iq2Inv);

			auto states2 = innerBandStructure.getState(q2);

			auto state2Energies = states2.getEnergies();
			auto nb2 = state2Energies.size();

			// if the meshes are the same (and gamma centered)
			// q3 will fall into the same grid, and it's easy to get
			if ( dontComputeQ3 ) {
				//	auto q3Plus = q1 + q2;
				//	auto q3Mins = q1 - q2; // this should be the one to use
				auto q3Plus = q1 + q2;
				auto q3Mins = q2 - q1;
				auto states3Plus = innerBandStructure.getState(q3Plus);
				auto states3Mins = innerBandStructure.getState(q3Mins);

				state3PlusEnergies = states3Plus.getEnergies();
				state3MinsEnergies = states3Mins.getEnergies();

				nb3Plus = state3PlusEnergies.size();
				nb3Mins = state3MinsEnergies.size();
				auto [cp,cm] = coupling3Ph->getCouplingSquared(
						states1, states2, states3Plus, states3Mins);
				couplingPlus = cp;
				couplingMins = cm;
//				couplingPlus = Eigen::Tensor<double,3>(nb1,nb2,nb3Plus);
//				couplingMins = Eigen::Tensor<double,3>(nb1,nb2,nb3Mins);
//				couplingPlus.setConstant(1.);
//				couplingMins.setConstant(1.);

				bose3PlusData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);
				bose3MinsData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);

				for ( long ib3=0; ib3<nb3Plus; ib3++ ) {
					auto ind3 = outerBandStructure.getIndex(q3Plus.getIndex(),
							ib3);
					bose3PlusData.col(ib3) = outerBose.data.col(ind3);
				}
				for ( long ib3=0; ib3<nb3Mins; ib3++ ) {
					auto ind3 = outerBandStructure.getIndex(q3Mins.getIndex(),
							ib3);
					bose3MinsData.col(ib3) = outerBose.data.col(ind3);
				}
			} else {
				// otherwise, q3 doesn't fall into the same grid
				// and we must therefore compute it from the hamiltonian

				Eigen::Vector3d q3PlusC = q1.getCoords(Points::cartesianCoords)
						+ q2.getCoords(Points::cartesianCoords);
				Eigen::Vector3d q3MinsC = q1.getCoords(Points::cartesianCoords)
						- q2.getCoords(Points::cartesianCoords);

				auto [eigvals3Plus,eigvecs3Plus] = h0->diagonalizeFromCoords(
						q3PlusC);
				auto [eigvals3Mins,eigvecs3Mins] = h0->diagonalizeFromCoords(
						q3MinsC);

				nb3Plus = eigvals3Plus.size();
				nb3Mins = eigvals3Mins.size();

				DetachedState states3Plus(q3PlusC, eigvals3Plus, numAtoms,
						nb3Plus, eigvecs3Plus);
				DetachedState states3Mins(q3MinsC, eigvals3Mins, numAtoms,
						nb3Mins, eigvecs3Mins);

				auto [cp,cm] = coupling3Ph->getCouplingSquared(states1,
						states2, states3Plus, states3Mins);
				couplingPlus = cp;
				couplingMins = cm;

				state3PlusEnergies = states3Plus.getEnergies();
				state3MinsEnergies = states3Mins.getEnergies();

				bose3PlusData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);
				bose3MinsData = Eigen::MatrixXd::Zero(numCalcs, nb3Plus);

				for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
					double temperature = statisticsSweep.getCalcStatistics(
							iCalc).temperature;
					for ( long ib3=0; ib3<nb3Plus; ib3++ ) {
						bose3PlusData(iCalc,ib3) = statistics.getPopulation(
								state3PlusEnergies(ib3), temperature);
					}
					for ( long ib3=0; ib3<nb3Mins; ib3++ ) {
						bose3MinsData(iCalc,ib3) = statistics.getPopulation(
								state3MinsEnergies(ib3), temperature);
					}
				}
			}

			for ( long ib1=0; ib1<nb1; ib1++ ) {
				en1 = state1Energies(ib1);
				ind1 = outerBandStructure.getIndex(iq1,ib1);
				if ( en1 < energyCutoff ) {
					continue;
				}

				for ( long ib2=0; ib2<nb2; ib2++ ) {
					en2 = state2Energies(ib2);
					ind2 = innerBandStructure.getIndex(iq2,ib2);
					if ( en2 < energyCutoff ) {
						continue;
					}

					// split into two cases since there may be different bands
					for ( long ib3=0; ib3<nb3Plus; ib3++ ) {
						en3Plus = state3PlusEnergies(ib3);
						if ( en3Plus < energyCutoff ) {
							continue;
						}

						double enProd = en1 * en2 * en3Plus;

						switch ( smearing->getType() ) {
						case ( DeltaFunction::gaussian ):
							deltaPlus = smearing->getSmearing(
									en1 + en2 - en3Plus);
							break;
						case ( DeltaFunction::adaptiveGaussian ):
							v1 = states1.getVelocity(ib1);
							v2 = states2.getVelocity(ib2);
							deltaPlus = smearing->getSmearing(
									en1 + en2 - en3Plus, v1+v2);
							break;
						default:
							deltaPlus = smearing->getSmearing(
									en3Plus-en1, iq2, ib2);
							break;
						}

//						deltaPlus = 1.;
//						couplingPlus(ib1,ib2,ib3) = 1.;

						// loop on temperature
						for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

							bose1 = outerBose.data(iCalc,ind1);
							bose2 = innerBose.data(iCalc,ind2);
							bose3Plus = bose3PlusData(iCalc,ib3);

							//Calculate transition probability W+
							ratePlus = pi * 0.25
									* bose1 * bose2 * ( bose3Plus + 1. )
									* couplingPlus(ib1,ib2,ib3)
									* deltaPlus / innerNumPoints / enProd;

							// note: to increase performance, we are in fact
							// using

							switch ( switchCase ) {
							case (0):
								// case of matrix construction
								// we build the scattering matrix S
								matrix(ind1,ind2) += ratePlus;
								linewidth->data(iCalc,ind1) += ratePlus;
								break;
							case (1):
								// case of matrix-vector multiplication
								// we build the scattering matrix A = S*n(n+1)
								for ( long i : {0,1,2} ) {
									outPopulation->data(3*iCalc+i,ind1) +=
											ratePlus *
											inPopulation->data(3*iCalc+i,ind2);
									outPopulation->data(3*iCalc+i,ind1) +=
											0.5 * ratePlus *
											inPopulation->data(3*iCalc+i,ind1);
								}
								break;
							case (2):
								// case of linewidth construction
								// there's still a missing norm done later
								linewidth->data(iCalc,ind1) += ratePlus;
								break;
							}
						}
					}

					for ( long ib3=0; ib3<nb3Mins; ib3++ ) {
						en3Mins = state3MinsEnergies(ib3);
						if ( en3Mins < energyCutoff ) {
							continue;
						}

						double enProd = en1 * en2 * en3Mins;

						switch ( smearing->getType() ) {
						case ( DeltaFunction::gaussian ):
							deltaMins1 = smearing->getSmearing(
									en1 + en3Mins - en2 );
							deltaMins2 = smearing->getSmearing(
									en2 + en3Mins - en1 );
							break;
						case ( DeltaFunction::adaptiveGaussian ):
							v1 = states1.getVelocity(ib1);
							v2 = states2.getVelocity(ib2);
							deltaMins1 = smearing->getSmearing(
									en1 + en3Mins - en2, v1-v2);
							deltaMins2 = smearing->getSmearing(
									en2 + en3Mins - en1, v1-v2);
							break;
						default:
							deltaMins1 = smearing->getSmearing(
									en1 - en3Mins, iq2, ib2);
							deltaMins2 = smearing->getSmearing(
									en1 - en3Mins, iq2, ib2);
							break;
						}

//						deltaMins1 = 1.;
//						deltaMins2 = 1.;
//						couplingMins(ib1,ib2,ib3) = 1.;

						for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {

							bose1 = outerBose.data(iCalc,ind1);
							bose2 = innerBose.data(iCalc,ind2);
							bose3Mins = bose3MinsData(iCalc,ib3);

							//Calculatate transition probability W-
							rateMins1 = pi * 0.25
									* bose3Mins * bose1 * ( bose2 + 1. )
									* couplingMins(ib1,ib2,ib3)
									* deltaMins1 / innerNumPoints / enProd;
							rateMins2 = pi * 0.25
									* bose2 * bose3Mins * ( bose1 + 1. )
									* couplingMins(ib1,ib2,ib3)
									* deltaMins2 / innerNumPoints / enProd;

							switch ( switchCase ) {
							case (0):
								// case of matrix construction
								matrix(ind1,ind2) -= rateMins1 + rateMins2;
								linewidth->data(iCalc,ind1) +=
										0.5 * rateMins2;
								break;
							case (1):
								// case of matrix-vector multiplication
								for ( long i : {0,1,2} ) {
									outPopulation->data(3*iCalc+i,ind1) +=
											- (rateMins1+rateMins2) *
											inPopulation->data(3*iCalc+i,ind2);
									outPopulation->data(3*iCalc+i,ind1) +=
											0.5 * (rateMins1+rateMins2) *
											inPopulation->data(3*iCalc+i,ind1);
								}
								break;
							case (2):
								// case of linewidth construction
								linewidth->data(iCalc,ind1) +=
										0.5 * rateMins2;
								break;
							}

						}
					}
				}
			}
		}
	}
	loopPrint.close();

	for ( long is=0; is<6*numPoints; is++ ) {
	std::cout << linewidth->data(0,is) * ryToCmm1 << "  !!!\n";
	}

	if ( switchCase == 0 ) {
		long numStates = outerBose.data.cols();
		long iCalc = 0;

		// case of matrix construction
		// we rescale the matrix in order to have S symmetrized
		for ( long ind1=0; ind1<numStates; ind1++ ) {
			bose1 = outerBose.data(iCalc,ind1);
			linewidth->data(iCalc,ind1) /= bose1;
			for ( long ind2=0; ind2<numStates; ind2++ ) {
				bose2 = outerBose.data(iCalc,ind2);
				matrix(ind1,ind2) /= sqrt(bose1*bose2);
			}
		}
	} else if ( switchCase == 2) {
		long numStates = outerBose.data.cols();
		long numCalcs = outerBose.data.rows();

		// case of linewidth construction
		// there's still a missing norm done later
		for ( long ind1=0; ind1<numStates; ind1++ ) {
			for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
				bose1 = outerBose.data(iCalc,ind1);
				linewidth->data(iCalc,ind1) /= bose1;
			}
		}
	}

	// we place the linewidths back in the diagonal of the scattering matrix
	// this because we need an MPI_allreduce on linewidth
	if ( switchCase == 0 ) { // case of matrix construction
		long iCalc = 0;
		for ( long i=0; i<outerBandStructure.getNumStates(); i++ ) {
			matrix(i,i) = linewidth->data(iCalc,i);
		}
	}



}
