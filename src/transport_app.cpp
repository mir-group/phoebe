#include <iostream>
#include <string>
#include "transport_app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "exceptions.h"

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
	auto [crystal, phononH0] =
			qeParser.parsePhHarmonic(context.getPhD2FileName());
	phononH0.setAcousticSumRule(context.getSumRuleD2());

	// Test
//	Eigen::VectorXd energies(3*crystal.getNumAtoms());
//	Eigen::Tensor<std::complex<double>,3> eigenvectors(3,crystal.getNumAtoms(),crystal.getNumAtoms()*3);
//	Eigen::VectorXd q(3);
//	q << 0.,0.,0.;
//	phononH0.diagonalize(q, energies, eigenvectors);
//	std::cout << energies.transpose();


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
