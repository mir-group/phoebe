#include <algorithm>
#include <exceptions.h>
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

IO::IO(int argc, char* argv[]) {

	char * inputFileName_ = getCmdOption(argv, argv + argc, "-in");
	char * outputFileName_ = getCmdOption(argv, argv + argc, "-out");
	if ( inputFileName_ == nullptr ) {
		Error e("Must provide an input file name on command line" ,1);
	}

	inputFileName = inputFileName_;
	if ( outputFileName_ != nullptr ) {
		outputFileName = outputFileName_;
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

std::string IO::getInputFileName() {
	return inputFileName;
}



void IO::welcome() {
	std::string welcomeMsg =
			"\n"
			"  8888888b.  888                        888\n"
			"  888   Y88b 888                        888\n"
			"  888    888 888                        888\n"
			"  888   d88P 88888b.   .d88b.   .d88b.  88888b.   .d88b.  \n"
			"  8888888P'  888 '88b d88''88b d8P  Y8b 888 '88b d8P  Y8b \n"
			"  888        888  888 888  888 88888888 888  888 88888888 \n"
			"  888        888  888 Y88..88P Y8b.     888 d88P Y8b.     \n"
			"  888        888  888  'Y88P'   'Y8888  88888P'   'Y8888' \n"
			"\n";
	std::cout << welcomeMsg;
}

void IO::goodbye() {
	std::cout << "Exiting program" << std::endl;
}
