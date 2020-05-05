#include "drift.h"

BulkTDrift::BulkTDrift(Context & context,
			ActiveBandStructure & activeBandStructure) :
			VectorBTE(context, activeBandStructure) {

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	Statistics statistics = activeBandStructure.getStatistics();

	for ( long is=0; is<numStates; is++ ) {
		double energy = activeBandStructure.getEnergy(is);
		Eigen::Vector3d velocity = activeBandStructure.getGroupVelocity(is);

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
			ActiveBandStructure & activeBandStructure) :
			VectorBTE(context, activeBandStructure) {

	Eigen::VectorXd temperatures = context.getTemperatures();
	Eigen::VectorXd chemicalPotentials = context.getChemicalPotentials();
	Statistics statistics = activeBandStructure.getStatistics();

	for ( long is=0; is<numStates; is++ ) {
		double energy = activeBandStructure.getEnergy(is);
		Eigen::Vector3d velocity = activeBandStructure.getGroupVelocity(is);

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
