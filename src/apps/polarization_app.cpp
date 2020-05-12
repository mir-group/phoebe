#include <string>
#include "polarization_app.h"
#include "exceptions.h"
#include "constants.h"
#include "periodic_table.h"

// Compute the electronic polarization using the Berry connection
void ElectronPolarizationApp::run(Context & context) {
	std::cout << "Starting electron polarization calculation" << std::endl;

	// Read the necessary input files
	auto [crystal, h0] = parser.parseElHarmonicWannier(context);

	// first we make compute the band structure on the fine grid
	FullPoints points(crystal, context.getKMesh());
	bool withVelocities = false;
	bool withEigenvectors = false;
	FullBandStructure bandStructure = h0.populate(points, withVelocities,
			withEigenvectors);

	// now we build the Berry connection

	Eigen::Tensor<double,3> berryConnection(points.getNumPoints(),
			h0.getNumBands(), 3);
	berryConnection.setZero();

	for ( long ik=0; ik<points.getNumPoints(); ik++ ) {
		auto point = points.getPoint(ik);
		std::vector<Eigen::MatrixXcd> thisBerryConnection = h0.getBerryConnection(point);
		for ( long ib=0; ib<h0.getNumBands(); ib++ ) {
			for ( long i=0; i<3; i++ ) {
				std::complex x = thisBerryConnection[i](ib,ib);
				berryConnection(ik,ib,i) = x.real();
			}
		}
	}

	// now we can compute the polarization

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	long nTemp = temperatures.size();
	long nMu = chemicalPotentials.size();

	Statistics statistics = h0.getStatistics();
	Eigen::Tensor<std::complex<double>,3> polarization(nTemp, nMu, 3);
	polarization.setZero();
	double population;

	for ( long ik=0; ik<points.getNumPoints(); ik++ ) {
		auto point = bandStructure.getPoint(ik);
		auto state = bandStructure.getState(point);
		auto energies = state.getEnergies();

		for ( long ib=0; ib<h0.getNumBands(); ib++ ) {
			double energy = energies(ib);

			for ( long it=0; it<nTemp; it++ ) {
				double temp = temperatures(it);

				for ( long imu=0; imu<nMu; imu++ ) {
					double mu = chemicalPotentials(imu);

					population = statistics.getPopulation(energy, temp, mu);

					for ( long i=0; i<3; i++ ) {
						polarization(it,imu,i) -=
								population * berryConnection(ik,ib,i);
					}

				}
			}
		}
	}

	for ( long it=0; it<nTemp; it++ ) {
		for ( long imu=0; imu<nMu; imu++ ) {
			for ( long i=0; i<3; i++ ) {
				polarization(it,imu,i) /= points.getNumPoints()
						* crystal.getVolumeUnitCell();
			}
		}
	}

	// now we add the ionic term

	PeriodicTable periodicTable;
	Eigen::MatrixXd atomicPositions = crystal.getAtomicPositions();
	std::vector<std::string> atomicNames = crystal.getAtomicNames();
	long numAtoms = crystal.getNumAtoms();
	for ( long ia=0; ia<numAtoms; ia++ ) {
		Eigen::Vector3d position = atomicPositions.row(ia);
		double charge = double(periodicTable.getIonicCharge(atomicNames[ia]));
		for ( long it=0; it<nTemp; it++ ) {
			for ( long imu=0; imu<nMu; imu++ ) {
				for ( long i=0; i<3; i++ ) {
					polarization(it,imu,i) += charge * position(i)
							/ crystal.getVolumeUnitCell();
				}
			}
		}
	}

	// Save results to file
//	std::ofstream outfile("./polarization.dat");
//	outfile << "# Phonon density of states: frequency[Cmm1], Dos[1/Ry]\n";
//	for ( long i=0; i<numEnergies; i++ ) {
//		outfile << energies[i] * ryToCmm1 << "\t" << dos[i] << "\n";
//	}
	std::cout << "Electron polarization computed" << std::endl;
}

