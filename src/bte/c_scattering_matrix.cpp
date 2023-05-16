#include "c_scattering_matrix.h"
#include "constants.h"
#include "helper_3rd_state.h"
#include "io.h"
#include "mpiHelper.h"
#include <cmath>

CoupledScatteringMatrix::CoupledScatteringMatrix(Context &context_,
                                      StatisticsSweep &statisticsSweep_,
                                      BaseBandStructure &innerBandStructure_, // phonon
                                      BaseBandStructure &outerBandStructure_, // electron
                                      Interaction3Ph *coupling3Ph_,
                                      InteractionElPh *couplingElPh_,
                                      ElectronH0 *electronH0_, 
                                      PhononH0 *phononH0_)
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_, outerBandStructure_),
      coupling3Ph(coupling3Ph_), couplingElPh(couplingElPh_), 
      phononH0(phononH0_), electronH0(electronH0_)  {

    if(!innerBandStructure_.getParticle().isPhonon() && !outerBandStructure_.getParticle().isElectron()) {
      Error("Developer error: Tried to create CMatrix with bandstructures of wrong particle type!");
    }
    // TODO should we force scattering matrix in memory? 

}

// TODO here would we like to call teh constructors of the other Smatrix types as well, or is this good? 
 CoupledScatteringMatrix::CoupledScatteringMatrix(Context &context_, 
                                StatisticsSweep &statisticsSweep_,
                                BaseBandStructure &innerBandStructure_,
                                BaseBandStructure &outerBandStructure_,
                                Interaction3Ph *coupling3Ph_,
                                InteractionElPh *couplingElPh_, 
                                PhononH0 *phononH0_, ElectronH0 *electronH0_) {
    : ScatteringMatrix(context_, statisticsSweep_, innerBandStructure_, outerBandStructure_),
      coupling3Ph(coupling3Ph_), couplingElPh(couplingElPh_), phononH0(phononH0_), electronH0(electronH0_) {


  // TODO may need to shift something in the parent scattering matrix constructor
  // wrt the size allocated for this matrix, which should be (nphononstates, nelectron states)^2

  // TODO overall seems like we should allow for phonon and electron bandstructures to 
  // supply an optional shift on state indexing. Perhaps somewhere in their mapping functions? 
  // OHH what about just on the iBTE or state indices? 
  // couldl construct one band structure first, then add that bandstructure's nstates 
  // as a "shift" parameter

  // TODO at first, we probably need to block symmetry in the construction of this matrix. 
  // If that is horrible, we could resurrect the symmetry unfolding of the MT calculation. 

  // 3 cases:
  // theMatrix and linewidth is passed: we compute and store in memory the
  // scattering matrix and the diagonal
  // inPopulation+outPopulation is passed: we compute the action of the
  //       scattering matrix on the in vector, returning outVec = sMatrix*vector
  // only linewidth is passed: we compute only the linewidths

  // TODO there will be a linewidth construction case here that you should block
  int switchCase = 0;
  if (theMatrix.rows() != 0 && linewidth != nullptr && inPopulations.empty() && outPopulations.empty()) {
    switchCase = 0;
  } else if (theMatrix.rows() == 0 && linewidth == nullptr && !inPopulations.empty() && !outPopulations.empty()) {
    switchCase = 1;
  } else if (theMatrix.rows() == 0 && linewidth != nullptr && inPopulations.empty() && outPopulations.empty()) {
    switchCase = 2;
  } else {
    Error("Developer error: CMatrix builder found a non-supported case.");
  }

  if ((linewidth != nullptr) && (linewidth->dimensionality != 1)) {
    Error("Developer error: The linewidths shouldn't have dimensionality!");
  }

  // add in the different scattering contributions -------------------

  // precompute the fermi and bose factors
  Eigen::MatrixXd innerBose = precomputeOccupations(innerBandStructure); 
  Eigen::MatrixXd outerFermi = precomputeOccupations(outerBandStructure); 

  // generate the points on which these processes will be computed

  // TODO need to convert this to give me pairs in the different quadrants. 
  // probably have to overwrite this function with a special one 

  // may just be able to shift things around.
  std::vector<std::tuple<std::vector<int>, int>> qPairIterator =
                                        getIteratorWavevectorPairs(switchCase);

  // here we call the function to add ph-ph scattering
  addPhPhScattering(this, context, inPopulations, outPopulations, 
                                  switchCase, qPairIterator, 
                                  innerBose, innerBose,
                                  innerBandStructure, innerBandStructure); 

  // here we call the function to add el-ph scattering
  addElPhScattering(this, context, inPopulations, outPopulations, 
                                  switchCase, qPairIterator, 
                                  outerFermi, // it only uses one of occupation 
                                  innerBandStructure, outerBandStructure); 

  // Isotope scattering
  if (context.getWithIsotopeScattering()) {
    addIsotopeScattering(this, context, inPopulations, outPopulations, 
                            switchCase, qPairIterator, 
                            innerBose, innerBose, 
                            innerBandStructure, innerBandStructure); 
  }
  // Add boundary scattering
  // TODO there will be an issue here for the CMatrix -- 
  // we can't set particle to be phonon or electron... 
  // Call this twice for each section of the diagonal, 
  // in one case handing it the phonon bands, in the other the electron bands. 

  // TODO one issue that will definitely arise is the computation state indices, 
  // which happens inside this class automatically. 
  // Need to provide some... shift, to the state indices, or something like that
  // maybe an optional "state range" or "state shift"?
  if (!std::isnan(context.getBoundaryLength())) {
    if (context.getBoundaryLength() > 0.) {
      // phonon boundary scattering
      addBoundaryScattering(this, context, inPopulations, outPopulations, 
                            switchCase, innerBandStructure); 
      // electron boundary scattering
      addBoundaryScattering(this, context, inPopulations, outPopulations, 
                            switchCase, outerBandStructure); 
  } 

  // TODO add phel scattering 
  //if(!context.getElphFileName().empty()) {
      // this is a weird case because it requires another band structure object, 
      // which we don't have access to in the regular phonon class here. 
      //
      // we could make it so that the phononElectron scattering object... which maybe 
      // will also compute the drag terms, does all the creating of the band structure and whatnot that
      // is currently gumming up phonon_transport_app... 

    // outline: 
    // 1) compute a new electron bandstructure here
    // 2) call the phonon_electron scattering method supplying it and our phonon bandstructure as well as this matrix.
  //

  // MPI reduce the distributed data now that all the scattering is accounted for 
  if (switchCase == 1) {
    for (auto & outPopulation : outPopulations) {
      mpi->allReduceSum(&outPopulation.data);
    }
  } else {
    mpi->allReduceSum(&linewidth->data);
    if(outputUNTimes) {
      mpi->allReduceSum(&internalDiagonalUmklapp->data);
      mpi->allReduceSum(&internalDiagonalNormal->data);
    }
  }
  // I prefer to close loopPrint after the MPI barrier: all MPI are synced here
  loopPrint.close();

  // Average over degenerate eigenstates.
  // we turn it off for now and leave the code if needed in the future << what does this mean? 
  if (switchCase == 2) {
    degeneracyAveragingLinewidths(linewidth);
    if(outputUNTimes) {
      degeneracyAveragingLinewidths(internalDiagonalUmklapp.get());
      degeneracyAveragingLinewidths(internalDiagonalNormal.get());
    }
  }

  // some phonons like acoustic modes at the gamma, with omega = 0,
  // might have zero frequencies, and infinite populations. We set those
  // matrix elements to zero.
  if (switchCase == 0) {
    // case of matrix construction
    if (context.getUseSymmetries()) {
      for (auto iBte1 : excludeIndices) {
        linewidth->data.col(iBte1).setZero();
        for (auto iBte2 : excludeIndices) {
          for (int i : {0, 1, 2}) {
            for (int j : {0, 1, 2}) {
              BteIndex iBte1Idx(iBte1);
              BteIndex iBte2Idx(iBte2);
              CartIndex iCart1(i);
              CartIndex iCart2(j);
              int iMat1 = getSMatrixIndex(iBte1Idx, iCart1);
              int iMat2 = getSMatrixIndex(iBte2Idx, iCart2);
              theMatrix(iMat1, iMat2) = 0.;
            }
          }
        }
      }
    } else {
      for (auto iBte1 : excludeIndices) {
        linewidth->data.col(iBte1).setZero();
        for (auto iBte2 : excludeIndices) {
          theMatrix(iBte1, iBte2) = 0.;
        }
      }
    }
  } else if (switchCase == 1) {
    // case of matrix-vector multiplication
    for (auto iBte1 : excludeIndices) {
      for (auto & outPopulation : outPopulations) {
        outPopulation.data.col(iBte1).setZero();
      }
    }

  } else if (switchCase == 2) {
    // case of linewidth construction
    for (auto iBte1 : excludeIndices) {
      linewidth->data.col(iBte1).setZero();
      // TODO may need to toss U and N indices here?
    }
  }

  // we place the linewidths back in the diagonal of the scattering matrix
  // this because we may need an MPI_allReduce on the linewidths
  if (switchCase == 0) { // case of matrix construction
    int iCalc = 0;
    if (context.getUseSymmetries()) {
      // numStates is defined in scattering.cpp as # of irrStates
      // from the outer band structure
      for (int iBte = 0; iBte < numStates; iBte++) {
        BteIndex iBteIdx(iBte);
        // zero the diagonal of the matrix
        for (int i : {0, 1, 2}) {
          CartIndex iCart(i);
          int iMati = getSMatrixIndex(iBteIdx, iCart);
          for (int j : {0, 1, 2}) {
            CartIndex jCart(j);
            int iMatj = getSMatrixIndex(iBteIdx, jCart);
            theMatrix(iMati, iMatj) = 0.;
          }
          theMatrix(iMati, iMati) += linewidth->operator()(iCalc, 0, iBte);
        }
      }
    } else {
      for (int is = 0; is < numStates; is++) {
        theMatrix(is, is) = linewidth->operator()(iCalc, 0, is);
      }
    }
  }
}



