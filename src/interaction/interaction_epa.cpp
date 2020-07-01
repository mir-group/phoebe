#include "interaction_epa.h"

//default constructor
InteractionEpa::InteractionEpa(Eigen::VectorXd & bandExtrema_,
                               Eigen::VectorXd & binSize_,
                               Eigen::VectorXi & numBins_,
                               Eigen::VectorXd & phFreqAverage_,
                               Eigen::Tensor<double,4> & elPhMatAverage_) : bandExtrema(bandExtrema_), binSize(binSize_), numBins(numBins_), phFreqAverage(phFreqAverage_), elPhMatAverage(elPhMatAverage_) {
    
}

//copy constructor
InteractionEpa::InteractionEpa(const InteractionEpa & that) : bandExtrema(that.bandExtrema), binSize(that.binSize), numBins(that.numBins), phFreqAverage(that.phFreqAverage), elPhMatAverage(that.elPhMatAverage) {
    
}

//overload the assignment operator
InteractionEpa & InteractionEpa::operator=(const InteractionEpa & that) {
    
    //avoid self-assignment:
    if ( this != &that ) {
        bandExtrema = that.bandExtrema;
        binSize = that.binSize;
        numBins = that.numBins;
        phFreqAverage = that.phFreqAverage;
        elPhMatAverage = that.elPhMatAverage;
    }
    
    return *this;
}

Eigen::VectorXd InteractionEpa::getBandExtrema() {
    return bandExtrema;
}

Eigen::VectorXd InteractionEpa::getBinSize() {
    return binSize;
}

Eigen::VectorXi InteractionEpa::getNumBins() {
    return numBins;
}

Eigen::VectorXd InteractionEpa::getPhFreqAverage() {
    return phFreqAverage;
}

Eigen::Tensor<double,4> InteractionEpa::getElPhMatAverage() {
    return elPhMatAverage;
}

