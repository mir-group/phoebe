#include <string>
#include "phononH0.h"
#include "electron_h0_spline.h"

class QEParser {
public:
	std::tuple<Crystal, PhononH0> parsePhHarmonic(const std::string fileName);
	std::tuple<Crystal, FullPoints, ElectronH0Spline> parseElHarmonicSpline(
			const std::string fileName, double& fourierCutoff);
//	void parseElHarmonicWan(int argc, char** argv);
//	void parseElPhAnharmonicEPA(int argc, char** argv);
//	void parseElphAnharmonicWan(int argc, char** argv);
//	void parsePhAnharmonicD3(int argc, char** argv);
};
