#ifndef CONTEXT_H
#define CONTEXT_H

#include <limits>  // NaN
#include <string>
#include <vector>

#include "eigen.h"

/** Class containing the user input variables.
 * This class is mostly a container for the input variables.
 *
 * To add a new variable, write a get/set method, and modify setupFromInput().
 * The new variable should be set to
 */
class Context {
 private:
  std::string phD2FileName = "";
  std::string phD3FileName = "";
  std::string electronH0Name = "";
  std::string appName = "";
  std::string sumRuleD2 = "";
  int smearingMethod = -1;
  double smearingWidth = std::numeric_limits<double>::quiet_NaN();
  Eigen::VectorXd temperatures;
  std::vector<std::string> solverBTE;
  double convergenceThresholdBTE = 1e-5;
  long maxIterationsBTE = 50;

  bool scatteringMatrixInMemory = true;

  std::string windowType = "nothing";
  Eigen::Vector2d windowEnergyLimit;
  double windowPopulationLimit = std::numeric_limits<double>::quiet_NaN();

  Eigen::VectorXd dopings;
  Eigen::VectorXd chemicalPotentials = Eigen::VectorXd::Zero(1);
  double electronFourierCutoff = std::numeric_limits<double>::quiet_NaN();

  Eigen::Vector3i qMesh;
  Eigen::Vector3i kMesh;

  double fermiLevel = std::numeric_limits<double>::quiet_NaN();
  double numOccupiedStates = std::numeric_limits<double>::quiet_NaN();
  bool hasSpinOrbit = false;

  long dimensionality = 3;

  double dosMinEnergy = std::numeric_limits<double>::quiet_NaN();
  double dosMaxEnergy = std::numeric_limits<double>::quiet_NaN();
  double dosDeltaEnergy = std::numeric_limits<double>::quiet_NaN();

  Eigen::MatrixXd inputAtomicPositions;
  Eigen::VectorXi inputAtomicSpecies;
  std::vector<std::string> inputSpeciesNames;

  Eigen::Tensor<double, 3> pathExtrema;
  double deltaPath = 0.05;

  double constantRelaxationTime = std::numeric_limits<double>::quiet_NaN();
  bool withIsotopeScattering = true;  // add isotopes in phonon scatt matrix
  Eigen::VectorXd massVariance;       // mass variance for isotope scattering

  // add RTA boundary scattering in phonon scatt matrix
  // boundary length for isotope scattering
  double boundaryLength = std::numeric_limits<double>::quiet_NaN();

  // utilities for parsing

  std::vector<std::string> &split(const std::string &s, char delim,
                                  std::vector<std::string> &elems);
  std::vector<std::string> split(const std::string &s, char delim);
//  bool lineHasParameter(const std::string &line);
//  std::string parseParameterName(const std::string &line);

  //  Setter and getter for all the variables above
 public:
  /** gets the name of the file containing the lattice force constants.
   * For Quantum Espresso, this is the path to the file produced by q2r.
   * @return x: the file path.
   */
  std::string getPhD2FileName();

  std::string getPhD3FileName();

  /** gets the name of the file containing the electronic band structure.
   * For Quantum Espresso, this is the path to the XML file.
   * @return path: the file path.
   */
  std::string getElectronH0Name();

  /** gets the value of the cutoff to be used for the Fourier interpolation
   * of the band structure.
   * @return r: the cutoff value.
   */
  double getElectronFourierCutoff();

  /** gets the type of calculation to be run.
   * @return x: the name of the calculation, e.g. "electron-phonon" or
   * "phonon-phonon".
   */
  std::string getAppName();

  /** gets the sum rule to be imposed on the lattice force constants.
   * @return x: the name of the sum rule, i.e. "simple" or "crystal".
   */
  std::string getSumRuleD2();

  /** gets the mesh of points for harmonic phonon properties.
   * @return path: an array with 3 integers representing the q-point mesh.
   */
  Eigen::Vector3i getQMesh();

  /** gets the mesh of points for harmonic electronic properties.
   * @return path: an array with 3 integers representing the k-point mesh.
   */
  Eigen::Vector3i getKMesh();

  /** gets the Window type to be used to filter out states that don't
   * contribute to transport.
   * @return path: an array with 3 integers representing the k-point mesh.
   * @param windowType: a string, which can take values "none", "energy",
   * or "population"
   */
  std::string getWindowType();

  /** gets the values of energy limits to be used with a window on energies.
   * @return x: a vector of 2 doubles representing the minimum and maximum
   * values of energies that will be used
   */
  Eigen::Vector2d getWindowEnergyLimit();

  /** gets the value of population above which a state is considered active.
   * i.e. the state will be used if its occupation number deviates from 0 or
   * 1 by at least this amount.
   * @return x: the <double> value of the population threshold.
   */
  double getWindowPopulationLimit();

  /** gets the value of chemical potentials (in Rydbergs) to be used in the
   * calculation of transport properties
   * @return x: the vector of values for chemical potentials
   */
  Eigen::VectorXd getChemicalPotentials();

  /** gets the value of chemical potentials (in Rydbergs) to be used in the
   * calculation of transport properties
   * @return x: the vector of values for chemical potentials
   */
  Eigen::VectorXd getDopings();

  /** gets the value of temperatures (in Rydbergs) to be used in the
   * calculation of transport properties
   * @return x: the vector of values for temperatures
   */
  Eigen::VectorXd getTemperatures();

  std::vector<std::string> getSolverBTE();

  double getConvergenceThresholdBTE();

  long getMaxIterationsBTE();

  long getDimensionality();

  double getDosMinEnergy();

  double getDosMaxEnergy();

  double getDosDeltaEnergy();

  // Wannier90 output doesn't contain the crystal information.
  // the user must then supplement it in the input
  // at least, if there's no phonon run
  // we may change the behavior in the future, parsing another file
  Eigen::MatrixXd getInputAtomicPositions();

  Eigen::VectorXi getInputAtomicSpecies();

  std::vector<std::string> getInputSpeciesNames();

  Eigen::Tensor<double, 3> getPathExtrema();

  double getDeltaPath();

  double getFermiLevel();
  void setFermiLevel(const double &x);

  double getNumOccupiedStates();
  void setNumOccupiedStates(const double &x);

  bool getHasSpinOrbit();
  void setHasSpinOrbit(const bool &x);

  int getSmearingMethod();

  double getSmearingWidth();

  double getConstantRelaxationTime();

  bool getScatteringMatrixInMemory();

  bool getWithIsotopeScattering();

  Eigen::VectorXd getMassVariance();

  double getBoundaryLength();

  /** Reads the user-provided input file and saves the input parameters
   * @param fileName: path to the input file
   */
  void setupFromInput(std::string fileName);
};

#endif
