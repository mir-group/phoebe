#ifndef CONTEXT_H
#define CONTEXT_H

#include <string>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Core>

using namespace std;

class Context {
private:
	std::string phD2FileName = "";
	std::string phD3FileName = "";
	std::string scratchFolder = "./out";
	std::string electronH0Name = "";
	std::string elPhFileName = "";
	std::vector<std::string> calculation = {""};
	std::string sumRuleD2 = "";
	std::string smearingType = "";
	double smearingWidth = 0.;
	std::vector<double> temperatures = {0.};
	std::vector<double> isotopeCoupling = {0.};
	std::vector<std::string> solverBTE = {"SMA"};
//	std::vector<std::string> secondSoundWavevector = {0.,0.,0.};
	std::vector<double> surfaceScatteringSize = {0.};
	std::string surfaceScatteringDirection = "";
	std::string transportRegime = "";
	double variationalConvergenceThreshold = 1e-5;
	double variationalMaxSteps = 1e-5;
	std::vector<double> atomicMasses = {0.};
	std::string windowType;
	std::vector<double> windowLimits = {0.,0.};
	std::string elPhInterpolationType = "";
	std::vector<double> dopings = {0.};
	std::vector<double> chemicalPotentials = {0.};
	double relaxationTime = 0.;
	bool usePolarCorrection = true;
	bool useWignerCorrection = true;
	double bandGap = 0.;
	double electronFourierCutoff = 0.;

//  Setter and getter for all the variables above
public:
	/** sets the name of the file containing the lattice force constants.
	 * For Quantum Espresso, this is the path to the file produced by q2r.
	 * @param x: the file path.
	 */
	void setPhD2FileName(std::string x);
	/** gets the name of the file containing the lattice force constants.
	 * For Quantum Espresso, this is the path to the file produced by q2r.
	 * @return x: the file path.
	 */
	std::string getPhD2FileName();

//	void setPhD3FileName(std::string x);
//	std::string getPhD3FileName();

//	void setScratchFolder(std::string x);
//	std::string getScratchFolder();

	/** sets the name of the file containing the electronic band structure.
	 * For Quantum Espresso, this is the path to the XML file.
	 * @param x: the file path.
	 */
	void setElectronH0Name(std::string path);
	/** gets the name of the file containing the electronic band structure.
	 * For Quantum Espresso, this is the path to the XML file.
	 * @return path: the file path.
	 */
	std::string getElectronH0Name();

	/** sets the value of the cutoff to be used for the Fourier interpolation
	 * of the band structure.
	 * @param r: the cutoff value. r is the radius up to which all lattice
	 * vectors (in real space) are used for the interpolation.
	 */
	void setElectronFourierCutoff(double r);
	/** gets the value of the cutoff to be used for the Fourier interpolation
	 * of the band structure.
	 * @return r: the cutoff value.
	 */
	double& getElectronFourierCutoff();

//	void setElPhFileName(std::string x);
//	std::string getElPhFileName();

//	void setCalculation(std::string x);
//	std::string getCalculation();

	/** sets the sum rule to be imposed on the lattice force constants.
	 * @param x: the name of the sum rule, i.e. "simple" or "crystal".
	 * where simple is a rescaling of diagonal elements, and crystal finds the
	 * closest matrix compliant with the translational rules
	 */
	void setSumRuleD2(std::string x);
	/** gets the sum rule to be imposed on the lattice force constants.
	 * @return x: the name of the sum rule, i.e. "simple" or "crystal".
	 */
	std::string getSumRuleD2();

//	void setSmearingType(std::string x);
//	std::string getSmearingType();
//
//	void setSmearingWidth(double x);
//	double getSmearingWidth();
//
//	void setTemperatures(double* x);
//	std::vector<double> getTemperatures();
//
//	void setIsotopeCoupling(double* x);
//	std::vector<double> getIsotopeCoupling();
//
//	void setSolverBTE(std::string* x);
//	std::vector<std::string> getSolverBTE();
//
//	void setSecondSoundWavevector(std::string* x);
//	std::vector<double> getSecondSoundWavevector();
//
//	void setSurfaceScatteringSize(double* x);
//	std::vector<double> getSurfaceScatteringSize();
//
//	void setSurfaceScatteringDirection(std::string x);
//	std::string getSurfaceScatteringDirection ();
//
//	void setTransportRegime(std::string x);
//	std::string getTransportRegime();
//
//	void setVariationalConvergenceThreshold(double x);
//	double getVariationalConvergenceThreshold();
//
//	void setVariationalMaxSteps(int x);
//	int getVariationalMaxSteps();
//
//	void setAtomicMasses(double* x);
//	std::vector<double> getAtomicMasses();
//
//	void setWindowType(std::string x);
//	std::string getWindowType();
//
//	void setWindowLimits(double* x);
//	std::vector<double> getWindowLimits();
//
//	void setElPhInterpolationType(std::string x);
//	std::string getElPhInterpolationType();
//
//	void setDopings(double* x);
//	std::vector<double> getDopings();
//
//	void setChemicalPotentials(double* x);
//	std::vector<double> getChemicalPotentials();
//
//	void setRelaxationTime(bool x);
//	bool getRelaxationTime();
//
//	void setUsePolarCorrection(bool x);
//	bool getUsePolarCorrection();
//
//	void setUseWignerCorrection(bool x);
//	bool getUseWignerCorrection();
//
//	void setBandGap(double x);
//	double getBanGap();

	/** Reads the user-provided input file and saves the input parameters
	 * @param fileName: path to the input file
	 */
	void setupFromInput(string fileName);
};

#endif
