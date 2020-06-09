#include "drift.h"

BulkTDrift::BulkTDrift(StatisticsSweep & statisticsSweep_,
		FullBandStructure<FullPoints> & bandStructure_,
		const long & dimensionality_) :
		VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {

	Statistics statistics = bandStructure.getStatistics();
	for ( long is=0; is<numStates; is++ ) {
		double energy = bandStructure.getEnergy(is);
		Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

		for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
			auto [imu,it,idim] = loc2Glob(iCalc);
			double vel = velocity(idim);
			auto calcStat = statisticsSweep.getCalcStatistics(it,imu);
			auto chemicalPotential = calcStat.chemicalPotential;
			auto temperature = calcStat.temperature;

			double x = statistics.getDndt(energy, temperature,
					chemicalPotential) * vel;
			data(iCalc,is) = x;
		}
	}
}

BulkEDrift::BulkEDrift(StatisticsSweep & statisticsSweep_,
		FullBandStructure<FullPoints> & bandStructure_,
		const long & dimensionality_) :
				VectorBTE(statisticsSweep_, bandStructure_, dimensionality_) {

	Statistics statistics = bandStructure.getStatistics();
	for ( long is=0; is<numStates; is++ ) {
		double energy = bandStructure.getEnergy(is);
		Eigen::Vector3d velocity = bandStructure.getGroupVelocity(is);

		for ( long iCalc=0; iCalc<numCalcs; iCalc++ ) {
			auto [imu,it,idim] = loc2Glob(iCalc);
			double vel = velocity(idim);
			auto calcStat = statisticsSweep.getCalcStatistics(it,imu);
			auto chemicalPotential = calcStat.chemicalPotential;
			auto temperature = calcStat.temperature;

			double x = statistics.getDnde(energy, temperature,
					chemicalPotential) * vel;
			data(iCalc,is) = x;
		}
	}
}
