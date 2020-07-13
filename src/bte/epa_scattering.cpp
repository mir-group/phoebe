#include "epa_scattering.h"
#include "delta_function.h"
#include "context.h"
#include "io.h"
#include "interaction_epa.h"
#include "particle.h"
#include "constants.h"
#include "vector_bte.h"
#include "statistics_sweep.h"
#include "epa_parser.h"
#include "bandstructure.h"

EpaScattering::EpaScattering(Context & context_, StatisticsSweep & statisticsSweep_,
                             FullBandStructure & fullBandStructure_): context(context_), statisticsSweep(statisticsSweep_), fullBandStructure(fullBandStructure_) {
}

EpaScattering::~EpaScattering(){
}

EpaScattering::EpaScattering(const EpaScattering &that) : context(that.context), statisticsSweep(that.statisticsSweep), fullBandStructure(that.fullBandStructure) {
    
}

EpaScattering& EpaScattering::operator=(const EpaScattering &that) {
    if (this != &that) {
        context = that.context;
        statisticsSweep = that.statisticsSweep;
        fullBandStructure = that.fullBandStructure;
    }
    return *this;
}

BaseVectorBTE EpaScattering::setup(Context & context, StatisticsSweep & statisticsSweep, FullBandStructure & fullBandStructure) {
    
    long numStates = fullBandStructure.getNumStates();
    
    /*If constant relaxation time is specified in input, we don't need to calculate
    EPA lifetimes*/
    double constantRelaxationTime = context.getConstantRelaxationTime();
    if (constantRelaxationTime > 0.) {
        BaseVectorBTE crtRate(statisticsSweep, numStates, 1);
        crtRate.setConst(1. / constantRelaxationTime);
        return crtRate;
    }
    
    auto hasSpinOrbit = context.getHasSpinOrbit();
    int spinFactor = 2;
    if (hasSpinOrbit)
        spinFactor = 1;
    
    auto particle = fullBandStructure.getParticle();
    
    if (particle.isPhonon())
        Error e("Electronic bandstructure has to be provided");
    
    long numPoints = fullBandStructure.getNumPoints();
    long numCalcs = statisticsSweep.getNumCalcs();
    
    std::cout << std::endl;
    std::cout << "Calculate electronic density of states. " << std::endl;
    TetrahedronDeltaFunction tetrahedra(fullBandStructure);
    
    double fermiLevel = context.getFermiLevel();
    double energyRange = context.getEnergyRange();
    double minEnergy = fermiLevel - energyRange;
    double maxEnergy = fermiLevel + energyRange;
    double energyStep = context.getEnergyStep();

    //in principle, we should add 1 to account for ends of energy interval
    //i will not do that, because will work with the centers of energy steps
    long numEnergies = (long) (maxEnergy-minEnergy)/energyStep;
    
    //energies at the centers of energy steps
    Eigen::VectorXd energies(numEnergies);
    
    for ( long i=0; i != numEnergies; ++i ) {
        //add 0.5 to be in the middle of the energy step
        energies(i) = (i + 0.5) * energyStep + minEnergy;
    }
    
    Eigen::VectorXd dos(numEnergies); // DOS initialized to zero
    dos.setZero();
    
    //calculate the density of states at the energies in energies vector
    for ( long i=0; i != numEnergies; ++i ) {
        dos(i) += tetrahedra.getDOS(energies(i));
    }
    
    // Write DOS to the file
    std::ofstream outfile;
    outfile.open("./electron_dos.dat");
    outfile << "# Electronic density of states: energy[Ry], Dos[1/Ry]" << std::endl;
    
    for ( long i=0; i != numEnergies; ++i ) {
        outfile << energies(i) << "\t"
        << dos[i]/energyRyToEv << "\n";
    }
    outfile.close();
    std::cout << "DOS is calculated and written to electron_dos.dat" << std::endl;
    
    //get vector containing averaged phonon frequencies per mode
    InteractionEpa couplingEpa = EpaParser::parseAvCouplings(context);
    
    Eigen::VectorXd phFreqAverage = couplingEpa.getPhFreqAverage();
    int numPhFreq = phFreqAverage.size();
    
    //phJump - contains the values of phFreqAverage/energyStep
    //defines in which step of electron energy grid the electron energy will go
    //after phonon absorption/emission
    Eigen::VectorXd phJump(phFreqAverage.size());
    
    for ( auto i = 0; i != phFreqAverage.size(); ++i) {
        phJump(i) = phFreqAverage(i)/energyStep;
    }
    
    int numBandGroups = couplingEpa.getNumBandGroups();
    Eigen::VectorXd extrema = couplingEpa.getBandExtrema();
    Eigen::VectorXi numBins = couplingEpa.getNumBins();
    int numBinsMax = numBins.maxCoeff();
    Eigen::Tensor<double,4> elPhMatElements = couplingEpa.getElPhMatAverage();
    Eigen::VectorXd binSize = couplingEpa.getBinSize();
    
    BaseVectorBTE epaRate(statisticsSweep, numEnergies, 1);
    LoopPrint loopPrint("to calculate EPA scattering rates", "pairs of temperatures and chemical potentials" , numCalcs);
    
    //loop over temperatures and chemical potentials
    for ( long iCalc = 0; iCalc != numCalcs; ++iCalc ) {
        loopPrint.update();
        double temperature = statisticsSweep.getCalcStatistics(iCalc).temperature;
        double chemPotential = statisticsSweep.getCalcStatistics(iCalc).chemicalPotential;
        
        //loop over energies
        for ( long iEnergy = 0; iEnergy != numEnergies; ++iEnergy ) {
            
            double scatRateTemp = 0.0;
            
            //iTest: make sure that by phonon absorption or emission we will not go outside of
            //the considered energy range
            int iTest = (int) phJump(numPhFreq);
            if ( iEnergy >= iTest && iEnergy <= numEnergies-iTest) {
                
                //loop over phonon frequencies
                for ( int iPhFreq = 0; iPhFreq != numPhFreq; ++iPhFreq ) {
                    //Bose-Einstein distribution
                    double nBose = 1. / (exp(phFreqAverage(iPhFreq)/temperature)-1);
                    if ( nBose < 0. ) {
                        nBose = 0.;
                    }
                    /*Fermi-Dirac distributions for  the cases when phonon is absorbed or emitted*/
                    double nFermiAbsorption = particle.getPopulation(energies[iEnergy] + phFreqAverage(iPhFreq),temperature,chemPotential);
                    double nFermiEmission = particle.getPopulation(energies[iEnergy] - phFreqAverage(iPhFreq),temperature,chemPotential);
                    
                    int iJump = (int) phJump(iPhFreq);
                    
                    double iInterp = phJump(iPhFreq) - (double) iJump;
                    
                    double dosAbsorption = 0;
                    double dosEmission = 0;
                    
                    dosAbsorption = dos(iEnergy + iJump)*(1.0 - iInterp) + dos(iEnergy + iJump + 1)*iInterp;
                    
                    dosEmission = dos(iEnergy - iJump -1)*iInterp + dos(iEnergy - iJump)*(1.0-iInterp);
                    
                    int iBandGroup = 0;
                    
                    if ( energies(iEnergy) <= extrema.sum()/numBandGroups ) {
                        iBandGroup = 0;
                    } else {
                        iBandGroup = 1;
                    }
                    
                    double iBinPos = (energies(iEnergy) - extrema(iBandGroup))/binSize(iBandGroup);
                    
                    iBinPos = std::max(iBinPos,1.0e-12);
                    iBinPos = std::min(iBinPos,numBins(iBandGroup) - 1.0e-12);
                    int intBinPos = (int) iBinPos;
                    
                    double elPhAbsorption = 0;
                    double elPhEmission = 0;
                    Eigen::VectorXd elPhAvTemp(numBinsMax);
                    elPhAvTemp.setZero();
                    
                    double gAbsorption = 0;
                    double gEmission = 0;
                    
                    if ( numBins(iBandGroup)==1 ) {
                        gAbsorption = elPhMatElements(iPhFreq,0,0,iBandGroup);
                        gEmission = elPhMatElements(iPhFreq,0,0,iBandGroup);
                    } else {
                        for ( int i = 0; i != numBinsMax; ++i ) {
                            elPhAvTemp(i) = elPhMatElements(iPhFreq,i,intBinPos,iBandGroup);
                    }
                    
                    double iAbs = (energies(iEnergy) + phFreqAverage(iPhFreq) - extrema(iBandGroup)) / binSize(iBandGroup);
                    iAbs = std::max(iAbs,1.0e-12);
                    iAbs = std::min(iAbs,numBins(iBandGroup) - 1.0e-12);
                    int iAbsInt = (int) iAbs;
                        
                    gAbsorption = elPhAvTemp(iAbsInt);
                    
                    double iEmis = (energies[iEnergy] - phFreqAverage(iPhFreq) - extrema(iBandGroup)) / binSize(iBandGroup);
                    iEmis = std::max(iAbs,1.0e-12);
                    iEmis = std::min(iAbs,numBins(iBandGroup) - 1.0e-12);
                    int iEmisInt = (int) iEmisInt;
                        
                    gEmission = elPhAvTemp(iEmisInt);
                        
                    }
                    
                    scatRateTemp += gAbsorption * (nBose + nFermiAbsorption) * dosAbsorption + gEmission * (nBose + 1 -nFermiEmission) * dosEmission;
                    
                }
                epaRate.data(iCalc,iEnergy) = scatRateTemp*rydbergSi*twoPi/spinFactor/hBarSi; //in 1/second
            }
        }
    }
    loopPrint.close();
    return epaRate;
}


