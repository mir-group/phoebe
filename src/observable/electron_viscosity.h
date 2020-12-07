#ifndef ELECTRONVISCOSITY_H
#define ELECTRONVISCOSITY_H

#include "drift.h"
#include "el_scattering.h"
#include "observable.h"

/** Object to compute and store the electron viscosity.
 * TODO: test if electron viscosity works with this class as well.
 */
class ElectronViscosity : public Observable {
public:
  /** Object constructor. Simply stores references to the inputs
   * @param statisticsSweep: object with info on the temperature and chemical
   * potential loops
   * @param crystal: object with the crystal structure
   * @param BaseBandStructure: object with the quasiparticle energies and
   * velocities computed on a mesh of wavevectors.
   */
  ElectronViscosity(Context &context_, StatisticsSweep &statisticsSweep_,
                    Crystal &crystal_, BaseBandStructure &bandStructure_);

  /** copy constructor
   */
  ElectronViscosity(const ElectronViscosity &that);

  /** copy assignment
   */
  ElectronViscosity &operator=(const ElectronViscosity &that);

  /** Compute the viscosity within the relaxation time approximation
   * Stores it internally.
   * @param n: the relaxation times.
   */
  virtual void calcRTA(VectorBTE &tau);

  /** Prints the viscosity to screen for the user.
   */
  void print();

  /** Outputs the quantity to a json file.
   * @param outFileName: string representing the name of the json file
   */
  void outputToJSON(std::string outFileName);

  /** Computes the viscosity from the scattering matrix eigenvectors.
   * Stores it internally.
   * @param vector0: VectorBTE object with the energy-conserving eigenvector.
   * @param relTimes: the VectorBTE object with relaxon relaxation times.
   * @oaram sMatrix: the scattering matrix.
   * @oaram eigenvectors: the eigenvectors of the scattering matrix above.
   */
  void calcFromRelaxons(Eigen::VectorXd &eigenvalues,
                        ParallelMatrix<double> &eigenvectors);

protected:
  virtual int whichType();
  BaseBandStructure &bandStructure;
};

#endif
