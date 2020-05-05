#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "context.h"
#include "exceptions.h"
#include "constants.h"

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

std::vector<double> parseDoubleList(std::vector<std::string> lines,
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

long parseLong(std::vector<std::string> lines, std::string pattern) {
	long x = 0;
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

std::vector<long> parseLongList(std::vector<std::string> lines,
		std::string pattern) {
	std::vector<long> x;
	long xTemp;
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
				throw "Error in parseLongList";
			}
			if ( pos2 == std::string::npos ) {
				throw "Error in parseLongList";
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

std::vector<std::string> parseStringList(std::vector<std::string> lines,
		std::string pattern) {
	std::vector<string> x;
	std::string xTemp;
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

				xTemp = token; // convert to integer
				x.push_back(xTemp);

			    s.erase(0, pos1 + delimeter.length());
			}
//			Must not forget the last element in the list
			xTemp = s;
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

void Context::setupFromInput(std::string fileName) {
	std::vector<std::string> lines;
	std::string line;

// open input file and read content
    ifstream infile(fileName);
	while (std::getline(infile, line)) {
		lines.push_back(line);
	}
    infile.close();

	try {
		std::string tmp = parseString(lines, "phD2FileName");
		setPhD2FileName(tmp);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		setSumRuleD2(parseString(lines, "sumRuleD2"));
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::string tmp = parseString(lines, "electronH0Name");
		setElectronH0Name(tmp);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		double dval = parseDouble(lines, "electronFourierCutoff");
		setElectronFourierCutoff(dval);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<long> vecMesh = parseLongList(lines, "qMesh");
		Eigen::Vector3i qMesh_;
		qMesh_(0) = vecMesh[0];
		qMesh_(1) = vecMesh[1];
		qMesh_(2) = vecMesh[2];
		setQMesh(qMesh_);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<long> vecMesh = parseLongList(lines, "kMesh");
		Eigen::Vector3i kMesh_;
		kMesh_(0) = vecMesh[0];
		kMesh_(1) = vecMesh[1];
		kMesh_(2) = vecMesh[2];
		setKMesh(kMesh_);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::string x = parseString(lines, "windowType");
		setWindowType(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<double> x = parseDoubleList(lines, "windowEnergyLimit");
		Eigen::Vector2d x_;
		// we just make sure to order it
		x_ << std::min(x[0],x[1]), std::max(x[0],x[1]);
		setWindowEnergyLimit(x_);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		double x = parseDouble(lines, "populationLimit");
		setWindowPopulationLimit(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<double> x = parseDoubleList(lines, "chemicalPotentials");
		Eigen::VectorXd x_(x.size());
		for ( long i=0; i<x.size(); i++ ) {
			x_(i) = x_[i] / energyRyToEv;
		}
		setChemicalPotentials(x_);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<double> x = parseDoubleList(lines, "temperatures");
		Eigen::VectorXd x_(x.size());
		for ( long i=0; i<x.size(); i++ ) {
			x_(i) = x_[i] / temperatureAuToSi;
		}
		setTemperatures(x_);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::string x = parseString(lines, "app");
		setAppName(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		std::vector<std::string> x = parseStringList(lines, "solverBTE");
		setSolverBTE(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		double x = parseDouble(lines, "convergenceThresholdBTE");
		setConvergenceThresholdBTE(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		long x = parseLong(lines, "maxIterationsBTE");
		setMaxIterationsBTE(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

	try {
		long x = parseLong(lines, "dimensionality");
		setDimensionality(x);
	} catch (ParameterNotFound& e) {} // Do nothing!

};

void Context::setPhD2FileName(std::string x) {
	phD2FileName = x;
}

std::string Context::getPhD2FileName() {
	return phD2FileName;
}

void Context::setSumRuleD2(std::string x) {
	sumRuleD2 = x;
}

std::string Context::getSumRuleD2() {
	return sumRuleD2;
}

void Context::setElectronH0Name(std::string x) {
	electronH0Name = x;
}

std::string Context::getElectronH0Name() {
	if ( electronH0Name == "" ) {
		Error e("Electronic H0 filename not set", 1);
	}
	return electronH0Name;
}

void Context::setElectronFourierCutoff(const double & x) {
	electronFourierCutoff = x;
}

double& Context::getElectronFourierCutoff() {
	if ( electronFourierCutoff == 0. ) {
		Error e("Electronic Fourier Cutoff not set", 1);
	}
	return electronFourierCutoff;
}

void Context::setAppName(const std::string & x) {
	appName = x;
}

std::string Context::getAppName() {
	return appName;
}

void Context::setQMesh(const Eigen::Vector3i & x) {
	qMesh = x;
}

Eigen::Vector3i Context::getQMesh() {
	return qMesh;
}

void Context::setKMesh(const Eigen::Vector3i & x) {
	kMesh = x;
}

Eigen::Vector3i Context::getKMesh() {
	return kMesh;
}


void Context::setWindowType(std::string x) {
	windowType = x;
}

std::string Context::getWindowType() {
	return windowType;
}

void Context::setWindowEnergyLimit(Eigen::Vector2d x) {
	windowEnergyLimit = x;
}

Eigen::Vector2d Context::getWindowEnergyLimit() {
	return windowEnergyLimit;
}

void Context::setWindowPopulationLimit(double x) {
	windowPopulationLimit = x;
}

double Context::getWindowPopulationLimit() {
	return windowPopulationLimit;
}

void Context::setChemicalPotentials(Eigen::VectorXd x) {
	chemicalPotentials = x;
}

Eigen::VectorXd Context::getChemicalPotentials() {
	return chemicalPotentials;
}

void Context::setTemperatures(Eigen::VectorXd x) {
	temperatures = x;
}

Eigen::VectorXd Context::getTemperatures() {
	return temperatures;
}

void Context::setNumValenceElectrons(long x) {
	numValenceElectrons = x;
}

long Context::getNumValenceElectrons() {
	return numValenceElectrons;
}

void Context::setHomo(double x) {
	homo = x;
}

double Context::getHomo() {
	return homo;
}

void Context::setSolverBTE(std::vector<std::string> x) {
	solverBTE = x;
}

std::vector<std::string> Context::getSolverBTE() {
	return solverBTE;
}

void Context::setConvergenceThresholdBTE(double x) {
	convergenceThresholdBTE = x;
}

double Context::getConvergenceThresholdBTE() {
	return convergenceThresholdBTE;
}

void Context::setMaxIterationsBTE(long x) {
	maxIterationsBTE = x;
}

long Context::getMaxIterationsBTE() {
	return maxIterationsBTE;
}

void Context::setDimensionality(long x) {
	dimensionality = x;
}

long Context::getDimensionality() {
	return dimensionality;
}





