#include <iostream>
#include <fstream>
#include "transport_epa_app.h"
#include "qe_input_parser.h"
#include "context.h"
#include "constants.h"
#include "exceptions.h"
#include "io.h"
#include "eigen.h"
#include "points.h"
#include "epa_parser.h"
#include "interaction_epa.h"
#include "particle.h"
#include "bandstructure.h"
#include "statistics_sweep.h"
#include "epa_scattering.h"
#include "vector_bte.h"
#include "crystal.h"
#include <cmath>
#include <string>

void TransportEpaApp::run(Context & context) {
    
    std::cout << "Setup EPA transport calculation" << std::endl;

    double fermiLevel = context.getFermiLevel();
    if ( std::isnan(fermiLevel) ) {
	Error e("Fermi energy must be provided for EPA calculation");
}
    std::cout << "Fermi level: " << fermiLevel <<std::endl;
    
    
    std::cout << "Point 1 reached" << std::endl;
    // Read necessary input: xml file of QE.
    //name of xml file should be provided in the input
    //electronFourierCutoff should be provided in the input (should it be the same as encut in DFT?)
    //Crystal crystal(directUnitCell, atomicPositions, atomicSpecies,
    //speciesNames, speciesMasses, dimensionality);
    //ElectronH0Fourier electronH0(crystal, coarsePoints, coarseBandStructure,
    //fourierCutoff);
    
    auto [crystal, electronH0] = QEParser::parseElHarmonicFourier(context);
    
    std::cout << "Point 2 reached" << std::endl;
    
    //Read and setup k-point mesh for interpolating bandstructure
    FullPoints fullPoints(crystal, context.getKMesh());
    bool withVelocities = true;
    bool withEigenvectors = true;
    FullBandStructure bandStructure = electronH0.populate(fullPoints, withVelocities, withEigenvectors);
    
    std::cout << "Point 3 reached" << std::endl;
    
    auto a = bandStructure.getNumBands(); // number of bands
    auto b = bandStructure.getNumStates(); // numBands*numPoints
    auto c = bandStructure.getNumPoints(); //number of k-points in fine mesh
    std::cout << "number of bands: " << a << std::endl;
    std::cout << "number of states: " << b << std::endl;
    std::cout << "number of points: " << c << std::endl;
    
   
        
    // set the chemical potentials to zero, load temperatures
    StatisticsSweep statisticsSweep(context, & bandStructure);
    
    BaseVectorBTE scatteringRates = EpaScattering::setup(context, statisticsSweep, bandStructure);
    
    // build/initialize the scattering matrix and the smearing
//    PhScatteringMatrix scatteringMatrix(context, statisticsSweep,
                                  //      bandStructure, bandStructure, &coupling3Ph);
//    scatteringMatrix.setup();
    
    // solve the BTE at the relaxation time approximation level
    // we always do this, as it's the cheapest solver and is required to know
    // the diagonal for the exact method.
    
//    std::cout << "\n";
//    std::cout << std::string(80, '-') << "\n";
//    std::cout << "\n";
//    std::cout << "Solving BTE within the relaxation time approximation.\n";
    
    // compute the phonon populations in the relaxation time approximation.
    // Note: this is the total phonon population n (n != f(1+f) Delta n)
    
 //   auto dimensionality = context.getDimensionality();
 //   BulkTDrift drift(statisticsSweep, bandStructure, dimensionality);
 //   VectorBTE phononRelTimes = scatteringMatrix.getSingleModeTimes();
 //   VectorBTE popRTA = drift * phononRelTimes;
    
    // compute the thermal conductivity
   // PhononThermalConductivity phTCond(statisticsSweep, crystal, bandStructure);
 //   phTCond.calcFromPopulation(popRTA);
//    phTCond.print();
    
    
    
    
  //  InteractionEpa interactionEpa = EpaParser::parseAvCouplings(context);
    
  //  Eigen::VectorXd phFreqAverage = interactionEpa.getPhFreqAverage();
    
//    for (auto x : phFreqAverage) {
 //       std::cout << x << " ";
 //   }
    
//    std::cout << std::endl;
    
    //If spin-orbit coupling is included, each band has the weight of 1, but we have 2 times more
    //bands. If there is no spin-orbit coupling, the weight of each band is 2. No spin-orbit
    //coupling  by default
    
//    auto hasSpinOrbit = context.getHasSpinOrbit();
//    int spinOrbit = 2;
//    if (hasSpinOrbit)
 //       spinOrbit = 1;

 //   std::cout << "spinOrbit=" << spinOrbit << std::endl;

    //add CRT option?
    //scissor shift?
    
  //  auto [crystal, electronH0] = parser.parseElHarmonicFourier(context);
    
    std::cout << "Starting EPA transport calculation" << std::endl;


    //If spin-orbit coupling is included, each band has the weight of 1, but we have 2 times more
    //bands. If there is no spin-orbit coupling, the weight of each band is 2
    //if (hasSpinOrbit)
    //    spinOrbit = 1;
    //else
      //  spinOrbit = 2;
    
    //Read EPA input
    
//    int numBandGroups, numPhFreq;
 //   Eigen::VectorXd bandExtrema(numBandGroups), binSize(numBandGroups);
    
 //   EpaParser::parseAvCouplings(context);
    
    
    std::cout << "EPA transport calculation is finished" << std::endl;
    
}

//void TransportEpaApp::checkRequirements(Context & context) {
//    throwErrorIfUnset(context.getEpaEFileName(), "epaEFileName");
//    throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
//    throwErrorIfUnset(context.getQMesh(), "kMesh");
//    throwErrorIfUnset(context.getDosMinEnergy(), "dosMinEnergy");
//    throwErrorIfUnset(context.getDosMaxEnergy(), "dosMaxEnergy");
//    throwErrorIfUnset(context.getDosDeltaEnergy(), "dosDeltaEnergy");
//    throwErrorIfUnset(context.getElectronFourierCutoff(),
//                      "electronFourierCutoff");
//}
