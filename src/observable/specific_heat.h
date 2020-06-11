#ifndef SPECIFICHEAT_H
#define SPECIFICHEAT_H

#include "observable.h"

/** Class for computing and storing the specific heat of a crystal.
 */
class SpecificHeat : public Observable {
public:
	/** Constructor method
	 * @param statisticsSweep: a StatisticsSweep object containing information
	 * on the temperature loop
	 * @param crystal: a Crystal object. Mostly used for the volume
	 * @param bandStructure: bandStructure to use for computing specific heat
	 */
	SpecificHeat(StatisticsSweep & statisticsSweep_,
			Crystal & crystal_,
			FullBandStructure<FullPoints> & bandStructure_);

	/** Copy constructor
	 */
	SpecificHeat(const SpecificHeat & that);

	/** Copy assignment operator
	 */
	SpecificHeat & operator = (const SpecificHeat & that);

	/** Computes the specific heat at all requested temperatures.
	 */
	virtual void calc();

	/** Prints to screen the specific heat at all desired temperatures.
	 */
	void print();

	/** returns the specific heat computed at a specific chemical potential
	 * and temperature. The indices are controlled by the statisticsSweep
	 * object.
	 * @param ChemPotIndex: strong typed index of chemical potential.
	 * @param TempIndex: strong typed index of temperature.
	 */
	const double & get(const ChemPotIndex & imu, const TempIndex & it);

protected:
	virtual int whichType();
	FullBandStructure<FullPoints> & bandStructure;
};

#endif
