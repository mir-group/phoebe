#include <assert.h>
#include <string>
#include <iostream>

class Error {
public:
	Error(const std::string & errMessage, const int & errCode = 1);
};

class Warning {
public:
	Warning(const std::string & errMessage);
};
