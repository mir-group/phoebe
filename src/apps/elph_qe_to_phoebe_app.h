#ifndef EL_PH_BLOCH_TO_WAN_APP_H
#define EL_PH_BLOCH_TO_WAN_APP_H

#include "app.h"
#include "electron_h0_wannier.h"
#include "phonon_h0.h"

/** Main driver for the transport calculation
 */
class ElPhQeToPhoebeApp : public App {
public:
  /** Executes the program contained in this class.
   */
  void run(Context &context) override;

  /** Check that user provided all necessary input
   */
  void checkRequirements(Context &context) override;

  // these are public because need testing
  /** In this function, we read the header file produced by the patched QE
   * version. The header contains some global and/or electronic information
   * coming from Quantum-ESPRESSO.
   *
   * @param crystal: class with the crystal info,
   * @param phoebePrefixQE: prefix of file,
   * e.g. "si" corresponds to the file "si.phoebe.****.dat"
   * @return tuple: with elements:
   *   * Eigen::Vector3i: mesh of q-points
   *   * Eigen::Vector3i: mesh of k-points
   *   * Eigen::MatrixXd(3,nk): kPoint cartesian coordinates in Bohr^-1
   *   * Eigen::MatrixXd(3,nq): qPoint cartesian coordinates in Bohr^-1
   *   * Eigen::Matrix(nb,nk): electronic energies
   *   * int : number of irreducible q-points.
   *   * numQEBands: number of KS states used by QE in the ph.x run
   *   * numElectrons: number of electrons in the unit cell
   *     (KS states>= electrons)
   *   * numSpin: number of spin components.
   */
  static std::tuple<Eigen::Vector3i, Eigen::Vector3i, Eigen::MatrixXd,
                    Eigen::MatrixXd, Eigen::MatrixXd, int, int, int, int>
  readQEPhoebeHeader(Crystal &crystal, const std::string &phoebePrefixQE);

  /** This is the driver for converting QE data to Phoebe data for Wannier
   * interpolation of the el-ph coupling.
   *
   * @param context: class with user input values.
   * @param crystal: class describing crystal unit cell.
   * @param phononH0: class with the phonon harmonic Hamiltonian.
   * @param kPoints: coarse grid of k-points used by QE.
   * @param qPoints: coarse grid of q-points used by QE. Equal to kPoints.
   * @param numQEBands: number of Kohn-Sham states used by QE.
   * @param numModes: number of phonon modes = 3*numAtoms.
   * @param numIrrQPoints: number of irreducible q-points computed by QE.
   * @param numElectrons: number of electrons in the crystal unit cell.
   * @param numSpin: number of spin components used by QE.
   * @param energies: matrix with electronic energies.
   * @param kGridFull: list of k-point cartesian coordinates in Bohr^-1.
   * @param kMesh: mesh of k-points used by QE.
   * @param qMesh: mesh of q-points used by QE.
   * @param runTests: if true (*default False) runs some sanity checks.
   */
  void postProcessingWannier(
      Context &context, Crystal &crystal, PhononH0 &phononH0, Points &kPoints,
      Points &qPoints, int numQEBands, int numModes, int numIrrQPoints,
      int numElectrons, int numSpin, const Eigen::MatrixXd &energies,
      const Eigen::MatrixXd &kGridFull, const Eigen::Vector3i &kMesh,
      const Eigen::Vector3i &qMesh, bool runTests = false);

protected:
  /** Transform the electron-phonon coupling computed by QE
   * from the Bloch to the Wannier representation. Follows roughly the procedure
   * described in [PRB 76, 165108 (2007)].
   * This is the expensive function of this app.
   *
   * @param elBravaisVectors: list of Bravais lattice vectors for the electronic
   * Fourier transform
   * @param phBravaisVectors: list of Bravais lattice vectors for the phonon
   * Fourier transform
   * @param gFull: electron-phonon coupling in the Bloch space, with size
   * gFull(nBands, nBands, nModes, nKPoints, nQPoints), where
   * in the notation of [PRB 76, 165108 (2007)], is
   * g(i,j,nu,k,q) = g_{i,j,nu}(k,q)
   *
   * @param uMatrices: the matrices for the rotation to the maximally localized
   * basis set of Wannier functions. Has size U(nBands,nWannier,k). Note that
   * there can be more bands than Wannier functions due to the disentanglement
   * procedure of Wannier90.
   * @param phEigenvectors: tensor (nu,mu,q) with the phonon eigenvectors,
   * where the first index is the index on cartesian and atomic basis,
   * the second is the index on phonon modes, and the third is on q-points.
   * @param kPoints: FullPoints class of k-points.
   * @param qPoints: FullPoints class of q-points.
   * @param crystal: class describing the crystal unit cell.
   * @param phononH0: class with the Phonon harmonic Hamiltonian.
   * @return g: the electron-phonon coupling in the Wannier representation.
   *   in the notation of [PRB 76, 165108 (2007)], is
   *   g(n,m,mu,R_p,R_e) = g_{m,n,mu}(R_e,R_p).
   *   Note that the first index refers to what must be Fourier transformed
   *   to k, whereas PRB refers that to the band at k+q.
   */
  Eigen::Tensor<std::complex<double>, 5> static blochToWannier(
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::MatrixXd &phBravaisVectors,
      Eigen::Tensor<std::complex<double>, 5> &gFull,
      const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
      const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
      Points &kPoints, Points &qPoints, Crystal &crystal, PhononH0 &phononH0);

  Eigen::Tensor<std::complex<double>, 5> BlochToWannierEfficient(
      Context &context, const Eigen::MatrixXd &energies,
      const Eigen::MatrixXd &kGridFull, const int &numIrrQPoints,
      const int &numQEBands, const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::MatrixXd &phBravaisVectors,
      const Eigen::Tensor<std::complex<double>, 3> &uMatrices, Points &kPoints,
      Points &qPoints, Crystal &crystal, PhononH0 &phononH0);

  /** Returns the rotation that moves the wave-function from the (entangled)
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
  static Eigen::Tensor<std::complex<double>, 3>
  setupRotationMatrices(const std::string &wannierPrefix, Points &fullPoints);

  /** Reads the electron-phonon coupling computed by QE on a coarse grid.
   *
   * Note: numQEBands != numBands. The first is the number of bands computes
   * in QE, while the latter is the number of bands passed to Wannier90.
   *
   * Note: we pass kGridFull because we need to align the list of kPoints used
   * by QE with the "list" of kPoints of Phoebe.
   *
   * @param context: class with user input.
   * @param numModes: number of phonon modes.
   * @param numBands: number of bands used in QE
   * @param numWannier: number of Wannier functions
   * @param kPoints: FullPoints class with kPoints
   * @param qPoints: FullPoints class with qPoints
   * @param kGridFull: list of kPoints in cartesian coordinates in QE
   * @param numIrrQPoints: number of irreducible qPoints used by QE
   * @param numQEBands: number of KS bands used by QE.
   * @param energies: matrix with electronic energies
   * @return tuple with:
   *   * gFull: electron-phonon coupling in the Bloch space, with size
   *     gFull(nBands, nBands, nModes, nKPoints, nQPoints), where
   *     in the notation of [PRB 76, 165108 (2007)], is
   *     g(i,j,nu,k,q) = g_{i,j,nu}(k,q)
   *   * complex tensor (nu,mu,q) with the phonon eigenvectors,
   *     where the first index is the index on cartesian and atomic basis,
   *     the second is the index on phonon modes, and the third is on q-points.
   *   * MatrixXd(nb,nq) containing phonon energies
   */
  std::tuple<Eigen::Tensor<std::complex<double>, 5>,
             Eigen::Tensor<std::complex<double>, 3>, Eigen::MatrixXd>
  readGFromQEFile(Context &context, const int &numModes, const int &numBands,
                  const int &numWannier, Points &kPoints, Points &qPoints,
                  const Eigen::MatrixXd &kGridFull, const int &numIrrQPoints,
                  const int &numQEBands, const Eigen::MatrixXd &energies);

  std::tuple<Eigen::Tensor<std::complex<double>, 5>,
             Eigen::Tensor<std::complex<double>, 3>, Eigen::MatrixXd,
             Eigen::MatrixXd>
  readChunkGFromQE(const int &iqIrr, Context &context, Points &kPoints,
                   const int &numModes, const int &numQEBands,
                   const Eigen::VectorXi &ikMap);

  /** This method compares the energies computed by qe2wannier90
   * and the energies of quantum espresso pw.x, to compute the offset between
   * the two sets (wannier90 may skip some core states).
   */
  static int computeOffset(const Eigen::MatrixXd &energies,
                           const std::string &wannierPrefix);

  // valid choices of interpolation, either "wannier" or "epa"
  std::vector<std::string> choices;

  /** Main driver for generating the EPA coupling. UNTESTED!
   */
  void epaPostProcessing(Context &context, Eigen::MatrixXd &elEnergies,
                         Points &kPoints, Points &qPoints,
                         const int &numElectrons, const int &numSpin,
                         const int &numModes, const int &numIrrQPoints,
                         const int &numQEBands,
                         const Eigen::MatrixXd &kGridFull);

  /** Test that the transform Wannier To Bloch of the electronic Hamiltonian
   * returns the original results computed by QE/Wannier90 on a coarse grid.
   */
  static void testElectronicTransform(
      Points &kPoints, const std::string &wannierPrefix,
      const Eigen::MatrixXd &elBravaisVectors,
      const Eigen::Tensor<std::complex<double>, 3> &uMatrices,
      const Eigen::VectorXd &elDegeneracies, ElectronH0Wannier &electronH0);

  /** Test that the transform Wannier To Bloch of the phonon harmonic
   * Hamiltonian returns the original results computed by QE on a coarse grid.
   */
  static void testPhononTransform(
      Crystal &crystal, PhononH0 &phononH0, Points &qPoints,
      const Eigen::Tensor<std::complex<double>, 3> &phEigenvectors,
      const Eigen::MatrixXd &phBravaisVectors,
      const Eigen::VectorXd &phDegeneracies, const Eigen::MatrixXd &phEnergies);

  /** Test that the transform Wannier To Bloch of the electron-phonon coupling
   * returns the original results computed by QE on a coarse grid, and that
   * it is correctly computed tby the class InteractionElPh.
   */
  static void testBackTransform(Context &context, PhononH0 &phononH0,
                                Points &kPoints, Points &qPoints,
                                ElectronH0Wannier &electronH0, Crystal &crystal,
                                Eigen::Tensor<std::complex<double>, 5> &gFull);

  /** This function writes the el-ph coupling to file. There are three modes
   * of writing file.
   * 1) plain text file. The slowest option, but we allow the user to not use
   * the HDF5 file format.
   * 2) HDF5 v1 (default and preferred). The el-ph tensor is written to a single
   * dataspace. if compiled with parallel HDf5, the write operation is done
   * in parallel. Note that we do the write operation in bunches of siez ~2Gb.
   * This is done to avoid a limitation of the MPI library, which has problem
   * when offsets to the tensor become > max(INT_32).
   * 3) HDF5 v2. The el-ph tensor is split along the R_e index into different
   * datasets. This format has been introduced because we suspect the HDF5
   * library has limitations when the electron-phonon tensor is very large.
   * We suspect this happens when the tensor size > max(UINT_32).
   */
  static void writeWannierCoupling(
      Context &context, Eigen::Tensor<std::complex<double>, 5> &gWannier,
      const int &numFilledWannier, const int &numSpin, const int &numModes,
      const int &numWannier, const Eigen::VectorXd &phDegeneracies,
      const Eigen::VectorXd &elDegeneracies,
      const Eigen::MatrixXd &phBravaisVectors,
      const Eigen::MatrixXd &elBravaisVectors, const Eigen::Vector3i &qMesh,
      const Eigen::Vector3i &kMesh);
};

#endif
