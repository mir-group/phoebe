#include "qe_input_parser.h"
#include "constants.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>  // to declare istringstream
#include <algorithm> // to use .remove_if

struct FileFormatNotRecognized : public std::exception {
	const char * what () const throw ()
    {
    	return "Error reading the file input parameter";
    }
};

std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);

   if ( delimiter == ' ' ) {
	   for (std::string s; tokenStream >> s; ) {
		   tokens.push_back(s);
	   }
   } else {
	   while (std::getline(tokenStream, token, delimiter)) {
		   token.erase(std::remove_if(token.begin(), token.end(), ::isspace),
				   token.end());
		   tokens.push_back(token);
	   }
   }

   return tokens;
}

void QEParser::parsePhHarmonic(std::string fileName) {
//  Here we read the dynamical matrix of interatomic force constants
//	in real space.
//	Since the file is typically small, we don't worry about memory management

	std::string line;
    std::vector<std::string> lineSplit;

// open input file
    std::ifstream infile(fileName);

//    this would read all content
//	std::vector<std::string> lines;
//	while (std::getline(infile, line)) {
//		lines.push_back(line);
//	}

//  First line contains ibrav, celldm and other variables

    std::getline(infile, line);
    lineSplit = split(line, ' ');

    int numElements = std::stoi(lineSplit[0]);
    int numAtoms = std::stoi(lineSplit[1]);
    int ibrav = std::stoi(lineSplit[2]);

    std::vector<double> celldm = {0.,0.,0.,0.,0.,0.};
    celldm[0] = std::stod(lineSplit[3]);
    celldm[1] = std::stod(lineSplit[4]);
    celldm[2] = std::stod(lineSplit[5]);
    celldm[3] = std::stod(lineSplit[6]);
    celldm[4] = std::stod(lineSplit[7]);
    celldm[5] = std::stod(lineSplit[8]);

	std::vector<std::vector<double>> unitCell(3, std::vector<double> (3, 0.));
    if ( ibrav == 0 ) {
//    	In this case, unitCell is written in the file, in angstroms
    	for ( int i=0; i<3; i++ ) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				unitCell[i][j] = std::stod(lineSplit[j]) / distanceRyToAng;
			}
    	};
    };


//  Next, we read the atomic species
    std::vector<std::string> speciesNames;
    std::vector<double> speciesMasses;
    for ( int i=0; i<numElements; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, '\'');
		speciesNames.push_back(lineSplit[1]);
		speciesMasses.push_back(std::stod(lineSplit[2]) / massRyToAmu );
    };


    //  we read the atomic positions
    std::vector<std::vector<double>> atomicPositions(numAtoms,
    		std::vector<double> (3,0.));
    std::vector<int> atomicSpecies(numAtoms, 0);
    for ( int i=0; i<numAtoms; i++ ) {
		std::getline(infile, line);
		lineSplit = split(line, ' ');
		atomicSpecies[i] = std::stoi(lineSplit[1]);
		atomicPositions[i][0] = std::stod(lineSplit[2]);
		atomicPositions[i][1] = std::stod(lineSplit[3]);
		atomicPositions[i][2] = std::stod(lineSplit[4]);
    }

//  Read if hasDielectric
	std::getline(infile, line);
	line.erase(std::remove_if(line.begin(), line.end(), ::isspace), line.end());
	bool hasDielectric;
	if ( line == "T" ) {
		hasDielectric = true;
	} else {
		hasDielectric = false;
	}

//	if there are the dielectric info, we can read dielectric matrix
//	and the Born charges
	if ( hasDielectric ) {
		std::vector<std::vector<double>> dielectricMatrix(3, std::vector<double> (3, 0.));
	    for ( int i=0; i<3; i++) {
			std::getline(infile, line);
			lineSplit = split(line, ' ');
			for ( int j=0; j<3; j++) {
				dielectricMatrix[i][j] = std::stod(lineSplit[j]);
			}
	    }

		std::vector<std::vector<std::vector<double>>> bornCharges(numAtoms,
				std::vector<std::vector<double>> (3,
						std::vector<double>(3,0.)) );
	    for ( int iAtom=0; iAtom < numAtoms; iAtom++ ) {
	    	std::getline(infile, line);
	    	for ( int i=0; i<3; i++ ) {
		    	std::getline(infile, line);
				lineSplit = split(line, ' ');
				for ( int j=0; j<3; j++ ) {
					bornCharges[iAtom][i][j] = std::stod(lineSplit[j]);
				}
	    	}
	    }
	}






//    std::cout << numElements << " , " << numAtoms <<  "!!!!\n";

    infile.close();
	return;
};


