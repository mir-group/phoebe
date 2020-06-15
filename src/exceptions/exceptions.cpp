#include <assert.h>
#include <string>
#include <iostream>
#include "exceptions.h"

Error::Error(const std::string & errMessage, const int & errCode) {
	if ( errCode != 0 ) {
		std::cout << "Error!" << std::endl;
		std::cout << errMessage << std::endl;
		exit(errCode);
	}
}

Warning::Warning(const std::string & errMessage) {
	std::cout << "WARNING: " << errMessage << std::endl;
}

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw () {
		return "Error reading the file input parameter";
	}
};
