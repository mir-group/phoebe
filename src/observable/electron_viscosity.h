#ifndef ELECTRON_VISCOSITY_H
#define ELECTRON_VISCOSITY_H

#include "drift.h"
#include "observable.h"

/** Object to compute and store the electron viscosity.
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
  void outputToJSON(const std::string& outFileName);

  /** Computes the viscosity from the scattering matrix eigenvectors.
   * Stores it internally.
   * @param vector0: VectorBTE object with the energy-conserving eigenvector.
   * @param relTimes: the VectorBTE object with relaxon relaxation times.
   * @param sMatrix: the scattering matrix.
   * @param eigenvectors: the eigenvectors of the scattering matrix above.
   */
  void calcFromRelaxons(Eigen::VectorXd &eigenvalues,
                        ParallelMatrix<double> &eigenvectors);

  /** Helper function to print information about the scalar products with the
   * special eigenvectors.
   * @param eigenvectors: eigenvectors of the scattering matrix
   * @param numRelaxons: the number of relaxons which have been calculated
   */
  void relaxonEigenvectorsCheck(ParallelMatrix<double>& eigenvectors, int& numRelaxons);


protected:
  int whichType() override;
  BaseBandStructure &bandStructure;
  int alpha0 = -1; // the index of the energy eigenvector, to skip it
  int alpha_e = -1; // the index of the charge eigenvector, to skip it

};

#endif
