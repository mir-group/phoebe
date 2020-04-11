#include <assert.h>
#include <string>
#include <iostream>
#include "exceptions.h"

Error::Error(std::string errMessage, int errCode) {
	if ( errCode != 0 ) {
		std::cout << "Error!" << std::endl;
		std::cout << errMessage << std::endl;
		assert(errCode != 0);
	}
}

Warning::Warning(std::string errMessage) {
	std::cout << errMessage;
}

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw () {
		return "Error reading the file input parameter";
	}
};
