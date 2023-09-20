#ifndef PHONON_VISCOSITY_H
#define PHONON_VISCOSITY_H

#include "drift.h"
#include "observable.h"
#include "ph_scattering_matrix.h"

/** Object to compute and store the phonon viscosity.
 * TODO: test if electron viscosity works with this class as well.
 */
class PhononViscosity : public Observable {
public:
  /** Object constructor. Simply stores references to the inputs
   * @param statisticsSweep: object with info on the temperature and chemical
   * potential loops
   * @param crystal: object with the crystal structure
   * @param BaseBandStructure: object with the quasiparticle energies and
   * velocities computed on a mesh of wavevectors.
   */
  PhononViscosity(Context &context_, StatisticsSweep &statisticsSweep_,
                  Crystal &crystal_, BaseBandStructure &bandStructure_);


  /** Compute the viscosity within the relaxation time approximation
   * Stores it internally.
   * @param n: the absolute phonon population.
   */
  virtual void calcRTA(VectorBTE &n);

  /** Prints the viscosity to screen for the user.
   */
  void print();

  /** Outputs the quantity to a json file.
   * @param outFileName: string representing the name of the json file
   */
  void outputToJSON(const std::string& outFileName);

  /** Computes the viscosity from the scattering matrix eigenvectors.
   * Following Simoncelli PRX 2020.
   * @param eigenvalues: the VectorBTE object with relaxon eigenvalues.
   * @param eigenvectors: the eigenvectors of the scattering matrix above.
   */
   // previously used vector0: VectorBTE object with the energy-conserving eigenvector.
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
};

#endif
