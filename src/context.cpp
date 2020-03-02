#include "context.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>

//TODO: it would be nice to have a smoother handling of read errors, with
//some informations provided to the user.
// Also, we must provide default values if key is not found in input

struct ParameterNotFound : public std::exception {
	const char * what () const throw ()
    {
    	return "Input parameter not found";
    }
};

bool lineHasPattern(std::string line, std::string pattern) {
	bool bx = false;

	std::string sep;
	std::string s;
	std::string substr2;
	size_t pos;

//	If line is empty, there's nothing for sure. Return false
	if ( line.empty() ) {
		return bx;
	}

	sep = "=";
	pos = line.find(sep);
	if ( string::npos == pos ) {
		throw "Delimeter not found";
	}

	s = line.substr(0,pos);
	s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

	if ( s == pattern ) {
		bx = true;
	}
	return bx;
}

bool parseBool(std::vector<std::string> lines, std::string pattern) {
	bool x;
	bool found = false;

	for ( std::string line : lines) {
		if ( lineHasPattern(line, pattern) ) {
			found = true;
			std::string delimeter = "=";
			size_t pos = line.find(delimeter);
			std::string s = line.substr(pos+1);

		    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
			s.erase(std::remove_if(s.begin(), s.end(),
					[](char c) { return !isalpha(c); } ), s.end());

		    if ( s == "false" ) {
		    	x = false;
		    } else if ( s == "true" ) {
		    	x = true;
		    } else if ( s == "0" ) {
		    	x = false;
		    } else if ( s == "1" ) {
		    	x = true;
		    } else {
		    	throw "Couldn't fix boolean value";
		    }
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};

double parseDouble(std::vector<std::string> lines, std::string pattern) {
	double x = 0.;
	bool found = false;

	for ( std::string line : lines) {
		if ( lineHasPattern(line, pattern) ) {
			std::string delimeter = "=";
			size_t pos = line.find(delimeter);
			std::string value = line.substr(pos+1);
			x = std::stod(value); // convert to double
			found = true;
			break;
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};

std::vector<double>  parseDoubleList(std::vector<std::string> lines,
		std::string pattern) {
	std::vector<double> x;
	double xTemp;
	bool found = false;
	std::string token, s;

	for ( std::string line : lines) {
		if ( lineHasPattern(line, pattern) ) {

			std::string delimeter;
			size_t pos1;
			size_t pos2;

			delimeter = "[";
			pos1 = line.find_first_of(delimeter);
			delimeter = "]";
			pos2 = line.find_last_of(delimeter);

			if ( pos1 == std::string::npos ) {
				throw "Error in parseDoubleList";
			}
			if ( pos2 == std::string::npos ) {
				throw "Error in parseDoubleList";
			}

			s = line.substr(pos1+1,pos2-pos1-1);
			delimeter = ",";

			while ((pos1 = s.find(delimeter)) != std::string::npos) {
			    token = s.substr(0, pos1);

				xTemp = std::stod(token); // convert to integer
				x.push_back(xTemp);

			    s.erase(0, pos1 + delimeter.length());
			}
//			Must not forget the last element in the list
			xTemp = std::stod(s);
			x.push_back(xTemp);

			found = true;
			break;
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};

int parseInt(std::vector<std::string> lines, std::string pattern) {
	int x = 0;
	bool found = false;

	for ( std::string line : lines) {
		if ( lineHasPattern(line, pattern) ) {
			std::string delimeter = "=";
			size_t pos = line.find(delimeter);
			std::string value = line.substr(pos+1);
			x = std::stoi(value); // convert to integer
			found = true;
			break;
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};

std::vector<int> parseIntList(std::vector<std::string> lines,
		std::string pattern) {
	std::vector<int> x;
	int xTemp;
	bool found = false;
	std::string token, s;

	for ( std::string line : lines) {
		if ( lineHasPattern(line, pattern) ) {

			std::string delimeter;
			size_t pos1;
			size_t pos2;

			delimeter = "[";
			pos1 = line.find_first_of(delimeter);
			delimeter = "]";
			pos2 = line.find_last_of(delimeter);

			if ( pos1 == std::string::npos ) {
				throw "Error in parseIntList";
			}
			if ( pos2 == std::string::npos ) {
				throw "Error in parseIntList";
			}

			s = line.substr(pos1+1,pos2-pos1-1);
			delimeter = ",";

			while ((pos1 = s.find(delimeter)) != std::string::npos) {
			    token = s.substr(0, pos1);

				xTemp = std::stoi(token); // convert to integer
				x.push_back(xTemp);

			    s.erase(0, pos1 + delimeter.length());
			}
//			Must not forget the last element in the list
			xTemp = std::stoi(s);
			x.push_back(xTemp);

			found = true;
			break;
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};


std::string parseString(std::vector<std::string> lines, std::string pattern) {
	std::string x = "";
	bool found = false;

	for ( std::string line : lines) {

		if ( lineHasPattern(line, pattern) ) {
			std::string delimeter;
			size_t pos1;
			size_t pos2;

			delimeter = "'";
			pos1 = line.find_first_of(delimeter);
			pos2 = line.find_last_of(delimeter);

			if ( pos1 == std::string::npos ) {

				delimeter = "\"";
				pos1 = line.find_first_of(delimeter);
				pos2 = line.find_last_of(delimeter);
				if ( pos1 == std::string::npos ) {
					throw "Couldn't solve string parsing of "
					"pattern: " + pattern;
				}
			}

			if ( pos1 == pos2 ) {
				throw "Error parsing string from user input. "
				"Pattern " + pattern;
			}
			x = line.substr(pos1+1,pos2-pos1-1);
			found = true;
			break;
		}
	}
	if ( not found ) {
		throw ParameterNotFound();
	}
	return x;
};

void Context::setQCoarseMesh(int* x) {
	for (int i=0; i<3; i++) {
		qCoarseMesh[i] = x[i];
	}
};

std::vector<int> Context::getQCoarseMesh() {
	return qCoarseMesh;
};

void Context::setupFromInput(std::string fileName) {
	std::vector<std::string> lines;
	std::string line;

// open input file and read content
    ifstream infile(fileName);
	while (std::getline(infile, line)) {
		lines.push_back(line);
	}
    infile.close();

// parse each variable

//    double x;
//	std::vector<double> w;
//    int z;
//    std::string y;
//	std::vector<int> q;
//	bool t;

//	x = parseDouble(lines, "x");
//	std::cout << x;
//	std::cout << "\n";
//	y = parseString(lines, "y");
//	std::cout << y;
//	std::cout << "\n";
//	z = parseInt(lines, "z");
//	std::cout << z;
//	std::cout << "\n";
//	q = parseIntList(lines, "qCoarseMesh");
//	std::cout << q[0] << q[1] << q[2];
//	std::cout << "\n";
//	w = parseDoubleList(lines, "w");
//	std::cout << w[0] << " , " << w[1] << " , " << w[2];
//	std::cout << "\n";
//	t = parseBool(lines, "t");
//	std::cout << t;
//	std::cout << "\n";

//	std::cout << "Sto qua\n";

	try {
		qCoarseMesh = parseIntList(lines, "qCoarseMesh");
	}
	catch (ParameterNotFound& e) {} // Do nothing!

	std::cout << qCoarseMesh[0] << "\n";
};
