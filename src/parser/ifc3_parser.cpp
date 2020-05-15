#include <iostream>
#include <fstream>
#include "ifc3_parser.h"
#include "eigen.h"

Interaction3Ph IFC3Parser::parseFromShengBTE(Context & context,
		Crystal & crystal) {

	auto fileName = context.getPhD3FileName();
	long tripletCount; // triplet counter

	// Open IFC3 file
	std::ifstream infile(fileName);

	long numTriplets;
	// Number of triplets
	infile >> numTriplets;

	// Allocate readables
	ifc3Tensor = Eigen::Tensor<double,4>(numTriplets,3,3,3);
	cellPositions = Eigen::Tensor<double,3>(numTriplets,2,3);
	displacedAtoms = Eigen::Tensor<long,2>(numTriplets,3);

	for ( long i=0; i<numTriplets; i++ ) {// loop over all triplets
		// Triplet counter
		infile >> tripletCount;

		// Read position of 2nd cell
		infile >> cellPositions(i,0,0) >> cellPositions(i,0,1)
				>> cellPositions(i,0,2);

		// Read position of 3rd cell
		infile >> cellPositions(i,1,0) >> cellPositions(i,1,1)
				>> cellPositions(i,1,2);

		//Convert cell positions from Ang to Bohr
		for ( long a : {0,1,2} ) {
			for ( long b : {0,1,2} ) {
				cellPositions(i,a,b) /= distanceBohrToAng;
			}
		}

		// Read triplet atom indices
		infile >> displacedAtoms(i,0) >> displacedAtoms(i,1)
				>> displacedAtoms(i,2);
		for ( long a : {0,1,2} ) {
			displacedAtoms(i,a) = displacedAtoms(i,a) - 1;
			// go to zero based indexing
		}

		// Read the 3x3x3 force constants tensor
		for ( long a : {0,1,2} ) {
			for ( long b : {0,1,2} ) {
				for ( long c : {0,1,2} ) {
					long tmp[3]; // temporary int holder
					infile >> tmp[0] >> tmp[1] >> tmp[2] >>
							ifc3Tensor(i,a,b,c);
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
