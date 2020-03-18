#include <assert.h>
#include <string>
#include <iostream>

class Error {
public:
	Error(std::string errMessage, int errCode);
};

class Warning {
public:
	Warning(std::string errMessage);
};
