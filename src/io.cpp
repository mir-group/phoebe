#include "io.h"

IO::IO(std::string inputFileName_, std::string outputFileName_) {
	inputFileName = inputFileName_;
	outputFileName = outputFileName_;
	if ( outputFileName != "" ) {
		outputFile.open(outputFileName);
		std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
		std::cout.rdbuf(outputFile.rdbuf()); //redirect std::cout to out.txt!
	}
};

IO::~IO() {
	if ( outputFile.is_open() ) {
		outputFile.close();
	}
}
