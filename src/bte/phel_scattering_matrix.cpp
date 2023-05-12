#include "phel_scattering_matrix.h"
#include "constants.h"
#include "io.h"
#include "mpiHelper.h"
#include "periodic_table.h"

// TODO... love this get el and get ph bandstructure thing. 
// Could we perhaps extend this to all the scattering matrix objects? 
// We could also write things like... get phInner and getPhOuter 
// instead of inner and outer bandstructures, would be way clearer

// helper function to generate kpoint pairs for phel scattering
// This generates irr points for the ik1, iq3 indices
std::vector<std::tuple<int, std::vector<int>>> PhElScatteringMatrix::getIrrWavevectorPairs() {

  std::vector<std::tuple<int, std::vector<int>>> pairIterator;

  // here I parallelize over ik1 which is the outer loop on q-points
  std::vector<int> k1Iterator = getElBandStructure().parallelIrrPointsIterator();

  // I don't parallelize the inner band structure, the inner loop
  // populate vector with integers from 0 to numPoints-1
  std::vector<int> q3Iterator = getPhBandStructure().irrPointsIterator();

  for (size_t ik1 : k1Iterator) {
    auto t = std::make_tuple(int(ik1), q3Iterator);
    pairIterator.push_back(t);
  }
  return pairIterator;
}

PhElScatteringMatrix::PhElScatteringMatrix(Context &context_,
                                           StatisticsSweep &statisticsSweep_,
                                           BaseBandStructure &elBandStructure_,
                                           BaseBandStructure &phBandStructure_,
                                           InteractionElPhWan &couplingElPhWan_,
                                           ElectronH0Wannier &electronH0_)
    : ScatteringMatrix(context_, statisticsSweep_, elBandStructure_, phBandStructure_),
      couplingElPhWan(couplingElPhWan_), electronH0(electronH0_) {

  isMatrixOmega = true;
  // automatically false, as phel scattering is only the diagonal
  highMemory = false;
}

// In the phononElectron case, we only compute the diagonal of the
// scattering matrix. Therefore, we compute only the linewidths
void PhElScatteringMatrix::builder(VectorBTE *linewidth,
                                   std::vector<VectorBTE> &inPopulations,
                                   std::vector<VectorBTE> &outPopulations) {
  (void) inPopulations;
  (void) outPopulations;

  if (linewidth == nullptr) {
    Error("builderPhEl found a non-supported case");
  }
  if (linewidth->dimensionality != 1) {
    Error("Linewidths shouldn't have dimensionality");
  }

  // compute the phonon electron lifetimes 
  addPhElScattering(); 

  // TODO could we compute boundary or isotope scattering here? 
  // I think they are diagonal terms, so this would work. 
  // 
  
  mpi->allReduceSum(&linewidth->data);

  // Average over degenerate eigenstates.
  // we turn it off for now and leave the code if needed in the future
  degeneracyAveragingLinewidths(linewidth);

}
