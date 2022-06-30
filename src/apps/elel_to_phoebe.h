#ifndef EL_EL_BLOCH_TO_WAN_APP_H
#define EL_EL_BLOCH_TO_WAN_APP_H

#include "app.h"
#include "electron_h0_wannier.h"

/** Main driver for the transport calculation
 */
class ElElToPhoebeApp : public App {
public:
  /** Executes the program contained in this class.
   */
  void run(Context &context) override;

  /** Check that user provided all necessary input
   */
  void checkRequirements(Context &context) override;

  // these are public because need testing
//  /** In this function, we read the header file produced by the patched QE
//   * version. The header contains some global and/or electronic information
//   * coming from Quantum-ESPRESSO.
//   *
//   * @param crystal: class with the crystal info,
//   * @param phoebePrefixQE: prefix of file,
//   * e.g. "si" corresponds to the file "si.phoebe.****.dat"
//   * @return tuple: with elements:
//   *   * Eigen::Vector3i: mesh of k-points
//   *   * Eigen::MatrixXd(3,nk): kPoint cartesian coordinates in Bohr^-1
//   *   * Eigen::Matrix(nb,nk): electronic energies
//   *   * numQEBands: number of KS states used by QE in the ph.x run
//   *   * numElectrons: number of electrons in the unit cell
//   *     (KS states>= electrons)
//   */
//  static std::tuple<Eigen::Vector3i, Eigen::MatrixXd, Eigen::MatrixXd, int, int>
//  readHeader(Crystal &crystal, const std::string &phoebePrefixQE);

protected:
  /** This function writes the el-el coupling to file.
   */
  static void writeWannierCoupling(
      Context &context, Eigen::Tensor<std::complex<double>, 7> &gWannier,
//      const int &numElectrons,
      const Eigen::VectorXd &elDegeneracies,
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::Vector3i &kMesh);
};

#endif
