#include <iostream>
#include <string>
#include "transport_app.h"
#include "qe_input_parser.h"
#include "context.h"

using namespace std;

//TransportApp::TransportApp() {};

void TransportApp::setup(int argc, char** argv) {
	string inputFileName;
	string outputFileName;

	if ( argc < 3 ) {
		throw std::invalid_argument("Invalid input line arguments. "
				"Must provide input file and output file.");
	}

	inputFileName = argv[1];
	outputFileName = argv[2];

	std::cout << inputFileName;
	std::cout << "\n";
	std::cout << outputFileName;
	std::cout << "\n";

//	Read user input file

	Context context;
	context.setupFromInput(inputFileName);

	std::cout << context.qCoarseMesh[0];
	std::cout << "\n";

//	Read the necessary input files

	QEParser qeParser;
	qeParser.parsePhHarmonic("SnSe.fc");


//	PhH0 phH0;
//	phH0.readFile();
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
