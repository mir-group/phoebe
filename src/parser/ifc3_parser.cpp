#include <iostream>
#include <fstream>
#include "ifc3_parser.h"
#include "eigen.h"
#include "constants.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif 

Interaction3Ph IFC3Parser::parse(Context &context, Crystal &crystal) {
    auto fileName = context.getPhD3FileName();

    // Open IFC3 file
    std::ifstream infile(fileName);

    if (not infile.is_open()) {
        Error e("D3 file not found", 1);
    }

    // ShengBTE has a single integer on the first line
    // QE has more (9). We use this to differentiate the formats

    std::string line;
    std::getline(infile, line);
    std::istringstream iss(line);
    std::string item;
    int counter = 0;
    while (iss >> item) {
        counter++;
    }

    if (counter == 1) {
        return parseFromShengBTE(context, crystal);
    } else {
        return parseFromQE(context, crystal);
    }
}

Interaction3Ph IFC3Parser::parseFromShengBTE(Context &context,
        Crystal &crystal) {

    auto fileName = context.getPhD3FileName();

    // Open IFC3 file
    std::ifstream infile(fileName);
    std::string line;

    if (not infile.is_open()) {
        Error e("D3 file not found");
    }

    // Number of triplets
    std::getline(infile, line);
    long numTriplets = std::stoi(line);

    // Allocate readables
    Eigen::Tensor<double, 4> ifc3Tensor(3, 3, 3, numTriplets);
    ifc3Tensor.setZero();
    Eigen::Tensor<double, 3> cellPositions(numTriplets, 2, 3);
    cellPositions.setZero();
    Eigen::Tensor<long, 2> displacedAtoms(numTriplets, 3);
    displacedAtoms.setZero();

    for (long i = 0; i < numTriplets; i++) {	// loop over all triplets

        // empty line
        std::getline(infile, line);
        // line with a counter
        std::getline(infile, line);

        // Read position of 2nd cell
        std::getline(infile, line);
        std::istringstream iss(line);
        std::string item;
        int j = 0;
        while (iss >> item) {
            cellPositions(i, 0, j) = std::stod(item) / distanceBohrToAng;
            j++;
        }

        // Read position of 3rd cell
        std::getline(infile, line);
        std::istringstream iss2(line);
        j = 0;
        while (iss2 >> item) {
            cellPositions(i, 1, j) = std::stod(item) / distanceBohrToAng;
            j++;
        }

        // Read triplet atom indices
        std::getline(infile, line);
        std::istringstream iss3(line);
        j = 0;
        long i0;
        while (iss3 >> i0) {
            displacedAtoms(i, j) = i0 - 1;
            j++;
        }

        // Read the 3x3x3 force constants tensor
        long i1, i2, i3;
        double d4;
        double conversion = pow(distanceBohrToAng, 3) / energyRyToEv;
        for (long a : { 0, 1, 2 }) {
            for (long b : { 0, 1, 2 }) {
                for (long c : { 0, 1, 2 }) {
                    std::getline(infile, line);
                    std::istringstream iss4(line);
                    while (iss4 >> i1 >> i2 >> i3 >> d4) {
                        ifc3Tensor(c, b, a, i) = d4 * conversion;
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

Interaction3Ph IFC3Parser::parseFromQE(Context &context, Crystal &crystal) {
    auto fileName = context.getPhD3FileName();

    const double threshold = 1.0e-20;

    // Open IFC3 file
    std::ifstream infile(fileName);
    std::string line;

    if (not infile.is_open()) {
        Error e("D3 file not found");
    }

    long numAtoms = crystal.getNumAtoms();
    long numSpecies = crystal.getSpeciesMasses().size();

    // The first few lines contain info on the crystal, which we ignore
    std::getline(infile, line);
    for (int i = 0; i < numSpecies; i++) {
        std::getline(infile, line);
    }
    for (int i = 0; i < numAtoms; i++) {
        std::getline(infile, line);
    }
    // line with dielectric
    std::getline(infile, line);

    // read the mesh of triplets
    std::getline(infile, line);

    long numTriplets = 0;

    // in this first loop we count the number of triplets that have
    // non zero derivative. This allows us to decrease the size of the matrix
    for (long na1 = 0; na1 < numAtoms; na1++) {
        for (long na2 = 0; na2 < numAtoms; na2++) {
            for (long na3 = 0; na3 < numAtoms; na3++) {
                for (long j1 : { 0, 1, 2 }) {
                    for (long j2 : { 0, 1, 2 }) {
                        for (long j3 : { 0, 1, 2 }) {
                            // suppress compiler warnings
                            (void) j1;
                            (void) j2;
                            (void) j3;

                            // read 6 integer which represent cartesian and
                            // atomic basis indeces
                            std::getline(infile, line);

                            // read the # of elements in this block
                            std::getline(infile, line);
                            std::istringstream iss2(line);
                            long nR;
                            iss2 >> nR;

                            for (long i = 0; i < nR; i++) {

                                std::getline(infile, line);
                                std::istringstream iss3(line);
                                int i1, i2, i3, i4, i5, i6;
                                double x1;
                                iss3 >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> x1;

                                if (abs(x1) > threshold) {
                                    numTriplets++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    infile.close();
    infile.clear();

    // now that we know the number of triplets, we build the matrix

    // Allocate readables
    Eigen::Tensor<double, 4> ifc3Tensor(3, 3, 3, numTriplets);
    ifc3Tensor.setZero();
    Eigen::Tensor<double, 3> cellPositions(numTriplets, 2, 3);
    cellPositions.setZero();
    Eigen::Tensor<long, 2> displacedAtoms(numTriplets, 3);
    displacedAtoms.setZero();

    // Open IFC3 file
    infile.open(fileName);

    // The first four lines contain info on the crystal, which we ignore
    std::getline(infile, line);
    for (int i = 0; i < numSpecies; i++) {
        std::getline(infile, line);
    }
    for (int i = 0; i < numAtoms; i++) {
        std::getline(infile, line);
    }
    std::getline(infile, line);

    // read the mesh of triplets
    std::getline(infile, line);

    long it = 0;

    for (long na1 = 0; na1 < numAtoms; na1++) {
        for (long na2 = 0; na2 < numAtoms; na2++) {
            for (long na3 = 0; na3 < numAtoms; na3++) {
                for (long j1 : { 0, 1, 2 }) {
                    for (long j2 : { 0, 1, 2 }) {
                        for (long j3 : { 0, 1, 2 }) {

                            // read 6 integer which represent cartesian and
                            // atomic basis indeces
                            std::getline(infile, line);

                            // read the # of elements in this block
                            std::getline(infile, line);
                            std::istringstream iss2(line);
                            long nR;
                            iss2 >> nR;

                            for (long i = 0; i < nR; i++) {

                                std::getline(infile, line);
                                std::istringstream iss3(line);
                                int i1, i2, i3, i4, i5, i6;
                                double x1;
                                iss3 >> i1 >> i2 >> i3 >> i4 >> i5 >> i6 >> x1;

                                if (abs(x1) > threshold) {
                                    cellPositions(it, 0, 0) = double(i1);
                                    cellPositions(it, 0, 1) = double(i2);
                                    cellPositions(it, 0, 2) = double(i3);

                                    cellPositions(it, 1, 0) = double(i4);
                                    cellPositions(it, 1, 1) = double(i5);
                                    cellPositions(it, 1, 2) = double(i6);

                                    // this is an index over the atomic basis
                                    // for the three atoms in the triplet
                                    displacedAtoms(it, 0) = na1;
                                    displacedAtoms(it, 1) = na2;
                                    displacedAtoms(it, 2) = na3;

                                    ifc3Tensor(j3, j2, j1, it) = x1;
                                    it++;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    infile.close();

    Interaction3Ph interaction3Ph(crystal, numTriplets, ifc3Tensor,
            cellPositions, displacedAtoms);

    return interaction3Ph;
}
