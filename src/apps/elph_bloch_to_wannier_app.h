#ifndef ELPHBLOCHTOWANAPP_H
#define ELPHBLOCHTOWANAPP_H

#include "app.h"

/** Main driver for the transport calculation
 */
class ElPhBlochToWannierApp : public App {
public:
  void run(Context &context);
  void checkRequirements(Context &context);
protected:
  /** Computes the transform of electron-phonon coupling from Bloch
   * to Wannier representation.
   */
  Eigen::Tensor<std::complex<double>, 5> blochToWannier(
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::MatrixXd &phBravaisVectors,
      Eigen::Tensor<std::complex<double>, 5> &g_full,
      const Eigen::Tensor<std::complex<double>,3> &uMatrices,
      const Eigen::Tensor<std::complex<double>,3> &phEigenvectors,
      FullPoints &kPoints, FullPoints &qPoints);

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
  Eigen::Tensor<std::complex<double>,3> setupRotationMatrices(
      const std::string &wannierPrefix, FullPoints &fullPoints);

};

#endif
