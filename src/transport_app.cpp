#include <iostream>
#include <string>
#include "transport_app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "points.h"
#include "io.h"

char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

std::tuple<std::string, std::string> TransportApp::parseInputArgs(int argc,
		char* argv[]) {
	char * inputFileName = getCmdOption(argv, argv + argc, "-in");
	char * outputFileName = getCmdOption(argv, argv + argc, "-out");
	if ( not inputFileName ) {
		Error e("Must provide an input file name on command line" ,1);
	}
	if ( outputFileName == nullptr ) {
		Error e("Must provide an output file name on command line" ,1);
	}
	return {inputFileName, outputFileName};
}

TransportApp::TransportApp(int argc, char** argv) {

auto [inputFileName, outputFileName] = parseInputArgs(argc, argv);

	IO io(inputFileName, outputFileName);

	std::cout << "Reading from input file: " << inputFileName << endl;

	// Read user input file

	Context context;
	context.setupFromInput(inputFileName);

	// Read the necessary input files

	QEParser qeParser;
	auto [crystal, phononH0] =
			qeParser.parsePhHarmonic(context.getPhD2FileName());
//  TODO: we could also use this syntax, if we fancy it
//	PhH0 phH0;
//	phH0.readFile();
	phononH0.setAcousticSumRule(context.getSumRuleD2());

	Eigen::VectorXd q(3);
	q << 0.,0.,0.1;
	auto [energies, eigenvectors] = phononH0.diagonalize(q);
	std::cout << energies.transpose() * ryToCmm1 << std::endl;

	auto [crystalEl, coarseElPoints, electronH0Spline] =
			qeParser.parseElHarmonicSpline(context.getElectronH0Name(),
					context.getElectronFourierCutoff());
	std::cout << "siamo usciti\n";

	Eigen::Vector3i mesh;
	mesh << 4, 4, 4;
	Points ps(crystal, mesh);

//
//	//	TODO: 'elph' shoudn't be a string, we should use a dictionary
//	//	and store which are the allowed values of calculations.
//
//	if ( context.getCalculation() == 'elph' ) {
//
//		if ( context.getInterpolation() == 'spline' ) {
//			ElH0Spline elH0Spline;
//			elH0Spline.readFile();
//
//		} else if ( context.getInterpolation() == "wannier" ) {
////			ElH0Wannier elH0Wannier;
////			elH0Wannier.readFile();
//			stopWithError("Not implemented");
//
//		} else {
//			stopWithError("key not recognized");
//		}
//	}

};
