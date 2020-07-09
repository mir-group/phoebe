#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

#include "eigen.h"
#include "constants.h"
#include "epa_parser.h"
#include "exceptions.h"
#include "interaction_epa.h"

InteractionEpa EpaParser::parseAvCouplings (Context & context) {
    //get the name of epa.e file
    auto fileName = context.getEpaEFileName();
    
    //open epa.e file for reading
    std::ifstream infile(fileName);
    
    //if infile cannot be opened, throw the error message and quit
    if ( !infile ) {
        Error e("epa.e file is not found", 1);
    }
    
    std::cout << "Reading " << fileName << " file" << std::endl;
    
    //Start reading infile
    std::string line;
    std::getline(infile, line);
    
    std::istringstream iss(line);
    
    //IMPORTANT: DOUBLECHECK HOW THIS WORKS IN CASE OF METALS
    //WE CAN TRY TO CHANGE THIS IMPLEMENTATION LATER ON WHEN EPA AVERAGING
    //WILL BE ADDED IN PHOEBE
    //numBandGroups: should be 2 (corresponds to valence and conduction bands)
    //numPhFreq: number of average phonon frequencies (equal to the number of phonon branches)
    int numBandGroups, numPhFreq;
    iss >> numBandGroups >> numPhFreq;
    
    //bandExtrema: vector of size 2:
    //bandExtrema(0) - top of the valence band (EV), bandExtrema(1) - bottom of the conduction band (EV)
    
    //binSize: vector of size 2:
    //binSize(0) - size of the energy bin for valence band (EV), binSize(1) - size of the energy bin for conduction band (EV)
    Eigen::VectorXd bandExtrema(numBandGroups), binSize(numBandGroups);
    bandExtrema.setZero();
    binSize.setZero();
    
    //numBins: vector of size 2:
    //numBins(0) - number of energy bins for valence band, numBins(1) - number of energy bins for conduction band
    Eigen::VectorXi numBins(numBandGroups);
    numBins.setZero();
    
    for (auto i = 0; i != numBandGroups; ++i) {
        
        std::getline(infile, line);
        std::istringstream iss1(line);
        
        iss1 >> bandExtrema(i) >> binSize(i) >> numBins(i);
    }
    
    //transform from Ev to Ry
    bandExtrema *= 1/energyRyToEv;
    binSize *= 1/energyRyToEv;
    
    //phFreqAverage - vector containing average phonon frequencies for numPhFreq branches (in cm^-1)
    Eigen::VectorXd phFreqAverage(numPhFreq);
    phFreqAverage.setZero();
    
    getline(infile,line);
    std::istringstream iss2(line);
    
    double temp;
    int i = 0;
    
    while ( iss2 >> temp ) {
        phFreqAverage(i) = temp;
        ++i;
    }
    
    //transform from Cm^-1 to Ry
    phFreqAverage /= ryToCmm1;

    auto numBinsMax = numBins.maxCoeff();

    //elPhMatAverage - tensor containing averaged squared electron-phonon matrix elements for each phonon branch
    //and each energy bin for valence and conduction bands
    Eigen::Tensor<double,4> elPhMatAverage(numPhFreq,numBinsMax,numBinsMax,numBandGroups);
    elPhMatAverage.setZero();
    
    for ( auto i = 0; i != numBandGroups; ++i ) {
        for ( auto j = 0; j != numBins(i); ++j ) {
            for ( auto k = 0; k != numBins(i); ++k ) {
               
                getline(infile,line);
                std::istringstream iss3(line);
                
                int ii;
                
                // we don't need the first 3 entries in each line, read them in ii variable
                iss3 >> ii >> ii >> ii;
                
                double temp;
                int l = 0;
                
                while (iss3 >> temp) {
                    elPhMatAverage(l,k,j,i) = temp/pow(energyRyToEv,2); // in Ry^2
                    ++l;
                }
                
            }
        }
    }
        
    infile.close();
    
    InteractionEpa interactionEpa(numBandGroups,bandExtrema,binSize,numBins,phFreqAverage,elPhMatAverage);
    return interactionEpa;
    
    std::cout << "... done" << std::endl;
    
}



