#include <string>
#include "phononH0.h"

class QEParser {
public:
//	void parseElHarmonicSpline(int argc, char** argv);
//	void parseElHarmonicWan(int argc, char** argv);
	std::tuple<Crystal, PhononH0> parsePhHarmonic(std::string fileName);

//	void parseElPhAnharmonicEPA(int argc, char** argv);
//	void parseElphAnharmonicWan(int argc, char** argv);
//	void parsePhAnharmonicD3(int argc, char** argv);
};
