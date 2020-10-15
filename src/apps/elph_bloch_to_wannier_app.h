#ifndef ELPHBLOCHTOWANAPP_H
#define ELPHBLOCHTOWANAPP_H

#include "app.h"
#include "phonon_h0.h"
#include "electron_h0_wannier.h"

/** Main driver for the transport calculation
 */
class ElPhQeToPhoebeApp : public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);

protected:
  /** Computes the transform of electron-phonon coupling from Bloch
   * to Wannier representation.
   */
  Eigen::Tensor<std::complex<double>, 5>
  blochToWannier(const Eigen::MatrixXd &elBravaisVectors,
                 const Eigen::MatrixXd &phBravaisVectors,
                 Eigen::Tensor<std::complex<double>, 5> &g_full,
                 const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
                 const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
                 FullPoints &kPoints, FullPoints &qPoints,
                 Crystal &crystal, PhononH0 &phononH0);

  /** Returns the rotation that moves the wavefunction from the (entangled)
   * Bloch representation to the disentangled Wannier representation
   *
   * @param wannierPrefix: path with the prefix used for the Wannier90
   * calculation (e.g. "si" or "./si")
   * @param fullPoints: coarse grid of wavevectors. Must be aligned with the
   * grid used in Wannier90 or an error will be raised.
   * @return U: tensor of shape (numBands, numWannier, numPoints), where
   * the wavevector index is aligned with fullPoints. If no disentanglement
   * was used in Wannier90, numWannier = numBands, otherwise,
   * numBands > numWannier due to the disentanglement. See Wannier90 docs.
   */
  Eigen::Tensor<std::complex<double>, 3>
  setupRotationMatrices(const std::string &wannierPrefix,
                        FullPoints &fullPoints);

  std::tuple<Eigen::Tensor<std::complex<double>, 5>,
             Eigen::Tensor<std::complex<double>, 3>, Eigen::MatrixXd>
  readGFromQEFile(Context &context, const int &numModes, const int &numBands,
                  const int &numWannier, FullPoints &kPoints,
                  FullPoints &qPoints, const Eigen::MatrixXd &kgridFull,
                  const int &numIrrQPoints, const int &numQEBands,
                  const Eigen::MatrixXd &energies);

  std::tuple<Eigen::Vector3i, Eigen::Vector3i, Eigen::MatrixXd, Eigen::MatrixXd,
             Eigen::MatrixXd, int, int, int, int>
  readQEPhoebeHeader(Crystal &crystal, const std::string &phoebePrefixQE);

  /** This method compares the energies computed by qe2wannier90
   * and the energies of quantum espresso pw.x, to compute the offset between
   * the two sets (wannier90 may skip some core states).
   */
  int computeOffset(const Eigen::MatrixXd &energies,
                    const std::string &wannierPrefix);

  std::vector<std::string> choices;

  void epaPostProcessing(Context &context,
                         Eigen::Tensor<std::complex<double>, 5> gFull,
                         Eigen::MatrixXd &elEnergies,
                         Eigen::MatrixXd &phEnergies, FullPoints &kPoints,
                         FullPoints &qPoints, const int &numElectrons,
                         const int &numSpin);

  void testElectronicTransform(
      Points &kPoints, const std::string &wannierPrefix,
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
      const Eigen::VectorXd &elDegeneracies, ElectronH0Wannier &electronH0);

  void testPhononTransform(
      Crystal &crystal, PhononH0 &phononH0, Points &qPoints,
      const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
      const Eigen::MatrixXd &phBravaisVectors,
      const Eigen::VectorXd &phDegeneracies, const Eigen::MatrixXd &phEnergies);

  };

#endif
