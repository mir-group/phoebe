#ifndef ELSCATTERING_H
#define ELSCATTERING_H

#include "interaction_elph.h"
#include "electron_h0_wannier.h"
#include "phonon_h0.h"
#include "scattering.h"
#include "vector_bte.h"

class ElScatteringMatrix : public ScatteringMatrix {
 public:
  ElScatteringMatrix(Context &context_, StatisticsSweep &statisticsSweep_,
                     BaseBandStructure &innerBandStructure_,
                     BaseBandStructure &outerBandStructure_,
                     PhononH0 &h0,
                     InteractionElPhWan *couplingElPhWan_ = nullptr
                     );

  ElScatteringMatrix(const ElScatteringMatrix &that);

  ElScatteringMatrix &operator=(const ElScatteringMatrix &that);

 protected:
  InteractionElPhWan *couplingElPhWan;
  PhononH0 &h0;

  double boundaryLength;
  bool doBoundary;

  virtual void builder(VectorBTE *linewidth,
                       std::vector<VectorBTE> &inPopulations,
                       std::vector<VectorBTE> &outPopulations);
};

#endif
