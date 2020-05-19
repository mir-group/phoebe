#include "drift.h"

BulkTDrift::BulkTDrift(Context & context,
		FullBandStructure<FullPoints> & bandStructure) :
			VectorBTE(context, bandStructure) {

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	Statistics statistics = activeBandStructure.getStatistics();

	for ( long is=0; is<numStates; is++ ) {
		double energy = bandStructure.getEnergy(is);
		Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

		for ( long idim=0; idim<dimensionality; idim++ ) {
			double vel = velocity(idim);
			for ( long imu=0; imu<numChemPots; imu++ ) {
				double chemicalPotential = chemicalPotentials(imu);
				for ( long it=0; it<numTemps; it++ ) {
					double temperature = temperatures(it);
					long ic = glob2Loc(imu, it, idim);
					double x = statistics.getDndt(energy, temperature,
							chemicalPotential) * vel;
					data(ic,is) = x;
				}
			}
		}
	}
}

BulkEDrift::BulkEDrift(Context & context,
			FullBandStructure<FullPoints> & bandStructure) :
			VectorBTE(context, bandStructure) {

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	Statistics statistics = bandStructure.getStatistics();

	for ( long is=0; is<numStates; is++ ) {
		double energy = activeBandStructure.getEnergy(is);
		Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

		for ( long idim=0; idim<dimensionality; idim++ ) {
			double vel = velocity(idim);
			for ( long imu=0; imu<numChemPots; imu++ ) {
				double chemicalPotential = chemicalPotentials(imu);
				for ( long it=0; it<numTemps; it++ ) {
					double temperature = temperatures(it);
					long ic = glob2Loc(imu, it, idim);
					double x = statistics.getDnde(energy, temperature,
							chemicalPotential) * vel;
					data(ic,is) = x;
				}
			}
		}
	}
}
