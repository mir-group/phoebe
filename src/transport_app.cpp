#include <iostream>
#include <string>
#include "transport_app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "points.h"

using namespace std;

//TransportApp::TransportApp() {};

void TransportApp::setup(int argc, char** argv) {
	string inputFileName;
	string outputFileName;

	if ( argc < 5 ) {
		Error e("Invalid input line arguments. "
				"Syntax is phoebe -in input -out output.", 1);
	}

	inputFileName = argv[2];
	outputFileName = argv[4];

	std::cout << "Reading from input file: " << inputFileName << endl;

	// Read user input file

	Context context;
	context.setupFromInput(inputFileName);

	// Read the necessary input files

	QEParser qeParser;

	// structured binding
	auto [crystal, phononH0] =
			qeParser.parsePhHarmonic(context.getPhD2FileName());
//  TODO: we could also use this syntax, if we fancy it
//	PhH0 phH0;
//	phH0.readFile();

	phononH0.setAcousticSumRule(context.getSumRuleD2());

	Eigen::Vector3i mesh;
	mesh << 4, 4, 4;
	Points ps(crystal, mesh);

	// Test
	Eigen::VectorXd q(3);
	q << 0.,0.,0.;
	auto [energies, eigenvectors] = phononH0.diagonalize(q);
	std::cout << energies.transpose() * ryToCmm1 << std::endl;
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
