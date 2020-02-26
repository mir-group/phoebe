#include <iostream>
#include "transport_app.h"
#include "context.h"
#include <string>

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

	Context context;
	context.setupFromInput(inputFileName);

	std::cout << context.qCoarseMesh;
	std::cout << "\n";
};
