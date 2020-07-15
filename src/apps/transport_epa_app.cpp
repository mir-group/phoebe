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
#include "delta_function.h"
#include "statistics_sweep.h"
#include "epa_scattering.h"
#include "vector_bte.h"
#include "crystal.h"
#include "utilities.h"
#include "onsager.h"

void TransportEpaApp::run(Context & context) {

    double fermiLevel = context.getFermiLevel();
    if ( std::isnan(fermiLevel) ) {
	Error e("Fermi energy must be provided for EPA calculation");
}
    // Read necessary input: xml file of QE.
    //name of xml file should be provided in the input
    //electronFourierCutoff should be provided in the input (should it be the same as encut in DFT?)
    //Crystal crystal(directUnitCell, atomicPositions, atomicSpecies,
    //speciesNames, speciesMasses, dimensionality);
    //ElectronH0Fourier electronH0(crystal, coarsePoints, coarseBandStructure,
    //fourierCutoff);
    
    auto [crystal, electronH0] = QEParser::parseElHarmonicFourier(context);
    
    //Read and setup k-point mesh for interpolating bandstructure
    FullPoints fullPoints(crystal, context.getKMesh());
    bool withVelocities = true;
    bool withEigenvectors = true;
    
    //Fourier interpolation of the electronic band structure
    FullBandStructure bandStructure = electronH0.populate(fullPoints, withVelocities, withEigenvectors);
    
    long numBands = bandStructure.getNumBands(); // number of bands
    long numStates = bandStructure.getNumStates(); // numBands*numPoints
    long numPoints = bandStructure.getNumPoints(true); //number of k-points in fine mesh
    Particle particle = bandStructure.getParticle();
    
    //calculate electronic band structure and energy projected velocity tensor
    
    
    // set temperatures, chemical potentials and carrier concentrations
    StatisticsSweep statisticsSweep(context, & bandStructure);
    
    int dimensionality = crystal.getDimensionality();
    
    double energyRange = context.getEnergyRange();
    double minEnergy = fermiLevel - energyRange;
    double maxEnergy = fermiLevel + energyRange;
    double energyStep = context.getEnergyStep();
    
    //in principle, we should add 1 to account for ends of energy interval
    //i will not do that, because will work with the centers of energy steps
    long numEnergies = (long) (maxEnergy-minEnergy)/energyStep;
    std::cout << "Num energies: " << numEnergies << std::endl;
    
    //energies at the centers of energy steps
    Eigen::VectorXd energies(numEnergies);
    
    for ( long i=0; i != numEnergies; ++i ) {
        //add 0.5 to be in the middle of the energy step
        energies(i) = (i + 0.5) * energyStep + minEnergy;
    }
    
    //Calculate EPA scattering rates
    BaseVectorBTE scatteringRates = EpaScattering::setup(context, statisticsSweep, bandStructure, energies);
    
    Eigen::Tensor<double, 3>  energyProjVelocity(dimensionality, dimensionality, numEnergies);
    energyProjVelocity.setZero();
    
    TetrahedronDeltaFunction tetrahedra(bandStructure);
    
    std::cout << "Start calculating energy projected velocity tensor" << std::endl;
    for ( long iEnergy = 0; iEnergy != numEnergies; ++iEnergy ) {
        for ( long iState = 0; iState != numStates; ++iState ) {
                    
            auto [ikC, ibC] = bandStructure.getIndex(iState);
            int ik = ikC.get();
            int ib = ibC.get();
            
            Eigen::Vector3d velocity = bandStructure.getGroupVelocity(iState);
            double deltaFunction = tetrahedra.getSmearing(energies(iEnergy), ik, ib);
                    
            for ( int iBeta = 0; iBeta != dimensionality; ++iBeta ) {
                for ( int iAlpha = 0; iAlpha != dimensionality; ++iAlpha) {

                    energyProjVelocity(iAlpha,iBeta,iEnergy) += velocity(iAlpha)*velocity(iBeta)*deltaFunction/numPoints;
                }
            }
        }
    }
    
    std::cout << "\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << "\n";
    
    OnsagerCoefficients EpaData(statisticsSweep, crystal, bandStructure, context);
    
    EpaData.calcFromEPA(scatteringRates, energyProjVelocity, energies, energyStep, particle);
    
    EpaData.calcTransportCoefficients();
    
    
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
