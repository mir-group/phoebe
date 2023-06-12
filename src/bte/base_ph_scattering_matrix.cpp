#include "base_ph_scattering_matrix.h"

BasePhScatteringMatrix::BasePhScatteringMatrix(Context &context_,
                                   StatisticsSweep &statisticsSweep_,
                                   BaseBandStructure &innerBandStructure_,
                                   BaseBandStructure &outerBandStructure_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_, outerBandStructure_) {

}

/*BaseBandStructure BaseElScatteringMatrix::initialElBandStructure() {
  return innerBandStructure; 
}*/

/*BaseBandStructure BaseElScatteringMatrix::finalElBandStructure() {
  return outerBandStructure; 
}*/

