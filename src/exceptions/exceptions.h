#include <assert.h>
#include <string>
#include <iostream>

void error(std::string errMessage, int errCode) {
	if ( errCode != 0 ) {
		std::cout << errMessage << std::endl;
		assert(errCode != 0);
	}
}

void warning(std::string errMessage) {
	std::cout << errMessage;
}

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw () {
		return "Error reading the file input parameter";
	}
};
