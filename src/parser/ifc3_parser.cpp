#include <iostream>
#include <fstream>
#include "ifc3_parser.h"
#include "eigen.h"
#include "constants.h"

Interaction3Ph IFC3Parser::parseFromShengBTE(Context & context,
		Crystal & crystal) {

	auto fileName = context.getPhD3FileName();

	// Open IFC3 file
	std::ifstream infile(fileName);
	std::string line;

	// Number of triplets
	std::getline(infile, line);
	long numTriplets = std::stoi(line);

	// Allocate readables
	Eigen::Tensor<double,4> ifc3Tensor(numTriplets,3,3,3);
	ifc3Tensor.setZero();
	Eigen::Tensor<double,3> cellPositions(numTriplets,2,3);
	cellPositions.setZero();
	Eigen::Tensor<long,2> displacedAtoms(numTriplets,3);
	displacedAtoms.setZero();

	for ( long i=0; i<numTriplets; i++ ) {// loop over all triplets

		// empty line
		std::getline(infile, line);
		// line with a counter
		std::getline(infile, line);

		// Read position of 2nd cell
		std::getline(infile, line);
		std::istringstream iss(line);
		std::string item;
		int j = 0;
		while ( iss >> item ) {
			cellPositions(i,0,j) = std::stod(item) / distanceBohrToAng;
			j++;
		}

		// Read position of 3rd cell
		std::getline(infile, line);
		std::istringstream iss2(line);
		j = 0;
		while ( iss2 >> item ) {
			cellPositions(i,1,j) = std::stod(item) / distanceBohrToAng;
			j++;
		}

		// Read triplet atom indices
		std::getline(infile, line);
		std::istringstream iss3(line);
		j = 0;
		long i0;
		while ( iss3 >> i0 ) {
			displacedAtoms(i,j) = i0 - 1;
			j++;
		}

		// Read the 3x3x3 force constants tensor
		long i1, i2, i3;
		double d4;
		double conversion = pow(distanceBohrToAng,3) / energyRyToEv;
		for ( long a : {0,1,2} ) {
			for ( long b : {0,1,2} ) {
				for ( long c : {0,1,2} ) {
					std::getline(infile, line);
					std::istringstream iss4(line);
					while ( iss4 >> i1 >> i2 >> i3 >> d4 ) {
						ifc3Tensor(i,a,b,c) = d4 * conversion;
					}
				}
			}
		}
	}

	// Close IFC3 file
	infile.close();

	//TODO Round cellPositions to the nearest lattice vectors

	Interaction3Ph interaction3Ph(crystal, numTriplets, ifc3Tensor,
			cellPositions, displacedAtoms);

	return interaction3Ph;
}
