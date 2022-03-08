#include "interaction_3ph.h"
#include "mpiHelper.h"
#include "common_kokkos.h"

Interaction3Ph::Interaction3Ph(Crystal &crystal, Eigen::Tensor<double, 5> &D3,
                               Eigen::MatrixXd &cellPositions2,
                               Eigen::MatrixXd &cellPositions3,
                               Eigen::Tensor<double, 3> &weights2,
                               Eigen::Tensor<double, 3> &weights3)
    : crystal_(crystal) {

  numAtoms = crystal_.getNumAtoms();
  numBands = numAtoms * 3;
  nr2 = cellPositions2.cols();
  nr3 = cellPositions3.cols();
  // Copy everything to kokkos views
  Kokkos::realloc(cellPositions2_k, nr2, 3);
  Kokkos::realloc(cellPositions3_k, nr3, 3);
  Kokkos::realloc(weights2_k, nr2, numAtoms, numAtoms);
  Kokkos::realloc(weights3_k, nr2, numAtoms, numAtoms);

  auto cellPositions2_h = Kokkos::create_mirror_view(cellPositions2_k);
  auto cellPositions3_h = Kokkos::create_mirror_view(cellPositions3_k);
  auto weights2_h = Kokkos::create_mirror_view(weights2_k);
  auto weights3_h = Kokkos::create_mirror_view(weights3_k);
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < nr2; i++) {
      cellPositions2_h(i, j) = cellPositions2(j, i);
    }
    for (int i = 0; i < nr3; i++) {
      cellPositions3_h(i, j) = cellPositions3(j, i);
    }
  }
  for (int j = 0; j < numAtoms; j++) {
    for (int k = 0; k < numAtoms; k++) {
      for (int i = 0; i < nr2; i++) {
        weights2_h(i, j, k) = weights2(i, j, k);
      }
      for (int i = 0; i < nr3; i++) {
        weights3_h(i, j, k) = weights3(i, j, k);
      }
    }
  }
  Kokkos::deep_copy(cellPositions2_k, cellPositions2_h);
  Kokkos::deep_copy(cellPositions3_k, cellPositions3_h);
  Kokkos::deep_copy(weights2_k, weights2_h);
  Kokkos::deep_copy(weights3_k, weights3_h);

  Kokkos::realloc(D3_k, numBands, numBands, numBands, nr3, nr2);
  Kokkos::realloc(D3PlusCached_k, numBands, numBands, numBands, nr3);
  Kokkos::realloc(D3MinsCached_k, numBands, numBands, numBands, nr3);
  auto D3_h = Kokkos::create_mirror_view(D3_k);
  for (int i1 = 0; i1 < numBands; i1++) {
    for (int i2 = 0; i2 < numBands; i2++) {
      for (int i3 = 0; i3 < numBands; i3++) {
        for (int i5 = 0; i5 < nr3; i5++) {
          for (int i4 = 0; i4 < nr2; i4++) {
            D3_h(i1, i2, i3, i5, i4) = D3(i1, i2, i3, i5, i4);
          }
        }
      }
    }
  }
  Kokkos::deep_copy(D3_k, D3_h);

  double memoryUsed = getDeviceMemoryUsage();
  kokkosDeviceMemory->addDeviceMemoryUsage(memoryUsed);
}

// copy constructor
Interaction3Ph::Interaction3Ph(const Interaction3Ph &that)
    : crystal_(that.crystal_), nr2(that.nr2), nr3(that.nr3),
      numAtoms(that.numAtoms), numBands(that.numBands) {
}

// assignment operator
Interaction3Ph &Interaction3Ph::operator=(const Interaction3Ph &that) {
  if (this != &that) {
    crystal_ = that.crystal_;
    nr2 = that.nr2;
    nr3 = that.nr3;
    numAtoms = that.numAtoms;
    numBands = that.numBands;
  }
  return *this;
}

Interaction3Ph::~Interaction3Ph() {
  double memoryUsed = getDeviceMemoryUsage(); // call this before deallocation
  kokkosDeviceMemory->removeDeviceMemoryUsage(memoryUsed);
  Kokkos::realloc(D3_k, 0, 0, 0, 0, 0);
  Kokkos::realloc(D3PlusCached_k, 0, 0, 0, 0);
  Kokkos::realloc(D3MinsCached_k, 0, 0, 0, 0);
  Kokkos::realloc(weights2_k, 0, 0, 0);
  Kokkos::realloc(weights3_k, 0, 0, 0);
  Kokkos::realloc(cellPositions2_k, 0, 0);
  Kokkos::realloc(cellPositions3_k, 0, 0);
}

void Interaction3Ph::cacheD3(const Eigen::Vector3d &q2_e) {
  // copy q2 to kokkos
  Kokkos::View<double *> q2("q2", 3);
  auto q2_h = Kokkos::create_mirror_view(q2);
  for (int i = 0; i < 3; i++) {
    q2_h(i) = q2_e(i);
  }
  Kokkos::deep_copy(q2, q2_h);

  Kokkos::complex<double> complexI(0.0, 1.0);

  // Need all variables to be local to be captured by lambda
  int nr2 = this->nr2;
  int nr3 = this->nr3;
  int numBands = this->numBands;

  auto D3PlusCached = this->D3PlusCached_k;
  auto D3MinsCached = this->D3MinsCached_k;
  auto cellPositions2 = this->cellPositions2_k;
  auto cellPositions3 = this->cellPositions3_k;
  auto weights2 = this->weights2_k;
  auto weights3 = this->weights3_k;
  auto D3 = this->D3_k;

  // precompute phases
  Kokkos::View<Kokkos::complex<double> *>
      phasePlus2("pp", nr2), phasePlus3("pp", nr3),
      phaseMins2("pp", nr2), phaseMins3("pp", nr3);

  Kokkos::parallel_for(
      "phase1loop", nr2, KOKKOS_LAMBDA(int ir2) {
        double argP = 0, argM = 0;
        for (int ic = 0; ic < 3; ic++) {
          argP += +q2(ic) * cellPositions2(ir2, ic);
          argM += -q2(ic) * cellPositions2(ir2, ic);
        }
        phasePlus2(ir2) = Kokkos::exp(complexI * argP);
        phaseMins2(ir2) = Kokkos::exp(complexI * argM);
      });

  Kokkos::parallel_for(
      "phase1loop", nr3, KOKKOS_LAMBDA(int ir3) {
        double argP = 0, argM = 0;
        for (int ic = 0; ic < 3; ic++) {
          argP += -q2(ic) * cellPositions3(ir3, ic);
          argM += +q2(ic) * cellPositions3(ir3, ic);
        }
        phasePlus3(ir3) = Kokkos::exp(complexI * argP);
        phaseMins3(ir3) = Kokkos::exp(complexI * argM);
      });

  // create cached D3
  Kokkos::parallel_for(
      "D3cacheloop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          {0, 0, 0, 0}, {numBands, numBands, numBands, nr3}),
      KOKKOS_LAMBDA(int ind1, int ind2, int ind3, int ir3) {
        // note: decompress2Indices doesn't compile like this on the GPU
        int at1 = ind1 / 3;
        int at2 = ind2 / 3;
        int at3 = ind3 / 3;
        // int at1 = std::get<0>(decompress2Indices(ind1,numAtoms,3));
        // int at2 = std::get<0>(decompress2Indices(ind2,numAtoms,3));
        // int at3 = std::get<0>(decompress2Indices(ind3,numAtoms,3));

        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int ir2 = 0; ir2 < nr2; ir2++) {// sum over all triplets
          tmpp += D3(ind1, ind2, ind3, ir3, ir2) * phasePlus2(ir2)
              * phasePlus3(ir3) * weights2(ir2, at1, at2) * weights3(ir3, at1, at3);
          tmpm += D3(ind1, ind2, ind3, ir3, ir2) * phaseMins2(ir2)
              * phaseMins3(ir3) * weights2(ir2, at1, at2) * weights3(ir3, at1, at3);
        }
        D3PlusCached(ind1, ind2, ind3, ir3) = tmpp;
        D3MinsCached(ind1, ind2, ind3, ir3) = tmpm;
      });
}

std::tuple<std::vector<Eigen::Tensor<double, 3>>,
           std::vector<Eigen::Tensor<double, 3>>>
Interaction3Ph::getCouplingsSquared(
    const std::vector<Eigen::Vector3d> &q1s_e, const Eigen::Vector3d &q2_e,
    const std::vector<Eigen::MatrixXcd> &ev1s_e, const Eigen::MatrixXcd &ev2_e,
    const std::vector<Eigen::MatrixXcd> &ev3Pluss_e,
    const std::vector<Eigen::MatrixXcd> &ev3Minss_e,
    const std::vector<int> &nb1s_e, const int nb2,
    const std::vector<int> &nb3Pluss_e, std::vector<int> &nb3Minss_e) {

  (void) q2_e;
  Kokkos::complex<double> complexI(0.0, 1.0);

  // Need all variables to be local to be captured by lambda
  int nr3 = this->nr3;
  int numBands = this->numBands;
  auto cellPositions2 = this->cellPositions2_k;
  auto cellPositions3 = this->cellPositions3_k;
  auto weights3 = this->weights3_k;
  auto D3 = this->D3_k;
  auto D3PlusCached = this->D3PlusCached_k;
  auto D3MinsCached = this->D3MinsCached_k;

  int nq1 = q1s_e.size();

  // MDRangePolicy loops are rectangular, need maximal dimensions
  int maxnb1 = *std::max_element(nb1s_e.begin(), nb1s_e.end()),
      maxnb3Plus = *std::max_element(nb3Pluss_e.begin(), nb3Pluss_e.end()),
      maxnb3Mins = *std::max_element(nb3Minss_e.begin(), nb3Minss_e.end());

  Kokkos::View<double **> q1s("q1s", nq1, 3);
  Kokkos::View<Kokkos::complex<double> **> ev2("ev2", nb2, numBands);
  Kokkos::View<Kokkos::complex<double> ***> ev1s("ev1s", nq1, maxnb1, numBands),
      ev3Pluss("ev3p", nq1, maxnb3Plus, numBands),
      ev3Minss("ev3m", nq1, maxnb3Mins, numBands);
  Kokkos::View<int *> nb1s("nb1s", nq1), nb3Pluss("nb3ps", nq1),
      nb3Minss("nb3ms", nq1);

  // copy everything to kokkos views
  {
    auto q1s_h = Kokkos::create_mirror_view(q1s);
    auto ev1s_h = Kokkos::create_mirror_view(ev1s);
    auto ev2_h = Kokkos::create_mirror_view(ev2);
    auto ev3Pluss_h = Kokkos::create_mirror_view(ev3Pluss);
    auto ev3Minss_h = Kokkos::create_mirror_view(ev3Minss);
    auto nb1s_h = Kokkos::create_mirror_view(nb1s);
    auto nb3Pluss_h = Kokkos::create_mirror_view(nb3Pluss);
    auto nb3Minss_h = Kokkos::create_mirror_view(nb3Minss);
    for (int i = 0; i < nq1; i++) {
      nb1s_h(i) = nb1s_e[i];
      nb3Pluss_h(i) = nb3Pluss_e[i];
      nb3Minss_h(i) = nb3Minss_e[i];
      for (int j = 0; j < 3; j++) {
        q1s_h(i, j) = q1s_e[i][j];
      }
      for (int j = 0; j < numBands; j++) {
        for (int k = 0; k < nb1s_e[i]; k++) {
          ev1s_h(i, k, j) = ev1s_e[i](j, k);
        }
      }
      for (int j = 0; j < nb3Pluss_e[i]; j++) {
        for (int k = 0; k < numBands; k++) {
          ev3Pluss_h(i, j, k) = ev3Pluss_e[i](k, j);
        }
      }
      for (int j = 0; j < nb3Minss_e[i]; j++) {
        for (int k = 0; k < numBands; k++) {
          ev3Minss_h(i, j, k) = ev3Minss_e[i](k, j);
        }
      }
    }
    for (int i = 0; i < numBands; i++) {
      for (int j = 0; j < nb2; j++) {
        ev2_h(j, i) = ev2_e(i, j);
      }
    }
    Kokkos::deep_copy(q1s, q1s_h);
    Kokkos::deep_copy(ev1s, ev1s_h);
    Kokkos::deep_copy(ev2, ev2_h);
    Kokkos::deep_copy(ev3Pluss, ev3Pluss_h);
    Kokkos::deep_copy(ev3Minss, ev3Minss_h);
    Kokkos::deep_copy(nb1s, nb1s_h);
    Kokkos::deep_copy(nb3Pluss, nb3Pluss_h);
    Kokkos::deep_copy(nb3Minss, nb3Minss_h);
  }

  Kokkos::View<Kokkos::complex<double> **> phases("pp", nq1, nr3);
  Kokkos::parallel_for(
      "tmpphaseloop",
      Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nq1, nr3}),
      KOKKOS_LAMBDA(int iq1, int ir3) {
        double arg = 0;
        for (int ic : {0, 1, 2}) {
          arg += -q1s(iq1, ic) * cellPositions3(ir3, ic);
        }
        phases(iq1, ir3) = exp(complexI * arg);
      });

  Kokkos::View<Kokkos::complex<double> ****> tmpPlus("tmpp", nq1, numBands,
                                                     numBands, numBands),
      tmpMins("tmpm", nq1, numBands, numBands, numBands);
  Kokkos::parallel_for(
      "tmploop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
          {0, 0, 0, 0}, {nq1, numBands, numBands, numBands}),
      KOKKOS_LAMBDA(int iq1, int iac1, int iac2, int iac3) {
        // note: decompress2Indices doesn't compile like this on the GPU
        int at1 = iac1 / 3;// std::get<0>(decompress2Indices(iac1,numAtoms,3));
        int at3 = iac3 / 3;// std::get<0>(decompress2Indices(iac3,numAtoms,3));
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int ir3 = 0; ir3 < nr3; ir3++) {// sum over all triplets
          tmpp += D3PlusCached(iac1, iac2, iac3, ir3) * phases(iq1, ir3) * weights3(ir3, at1, at3);
          tmpm += D3MinsCached(iac1, iac2, iac3, ir3) * phases(iq1, ir3) * weights3(ir3, at1, at3);
        }
        tmpPlus(iq1, iac1, iac2, iac3) = tmpp;
        tmpMins(iq1, iac1, iac2, iac3) = tmpm;
      });
  Kokkos::realloc(phases, 0, 0);

  Kokkos::View<Kokkos::complex<double> ****> tmp1Plus("t1p", nq1, maxnb1,
                                                      numBands, numBands),
      tmp1Mins("t1m", nq1, maxnb1, numBands, numBands);
  Kokkos::parallel_for(
      "tmp1loop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, numBands, numBands}),
      KOKKOS_LAMBDA(int iq1, int ib1, int iac2, int iac3) {
        int mask = ib1 < nb1s(iq1);
        Kokkos::complex<double> tmpp = 0, tmpm = 0;

        for (int iac1 = 0; iac1 < numBands; iac1++) {
          tmpp += tmpPlus(iq1, iac1, iac2, iac3) * ev1s(iq1, ib1, iac1);
          tmpm += tmpMins(iq1, iac1, iac2, iac3) * ev1s(iq1, ib1, iac1);
        }
        tmp1Plus(iq1, ib1, iac3, iac2) = tmpp * mask;
        tmp1Mins(iq1, ib1, iac3, iac2) = tmpm * mask;
      });
  Kokkos::realloc(tmpPlus, 0, 0, 0, 0);
  Kokkos::realloc(tmpMins, 0, 0, 0, 0);

  Kokkos::View<Kokkos::complex<double> ****> tmp2Plus("t2p", nq1, maxnb1, nb2,
                                                      numBands),
      tmp2Mins("t2m", nq1, maxnb1, nb2, numBands);
  Kokkos::parallel_for(
      "tmp2loop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, nb2, numBands}),
      KOKKOS_LAMBDA(int iq1, int ib1, int ib2, int iac3) {
        int mask = ib1 < nb1s(iq1);

        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int iac2 = 0; iac2 < numBands; iac2++) {
          tmpp += tmp1Plus(iq1, ib1, iac3, iac2) * ev2(ib2, iac2);
          tmpm += tmp1Mins(iq1, ib1, iac3, iac2) * Kokkos::conj(ev2(ib2, iac2));
        }
        tmp2Plus(iq1, ib1, ib2, iac3) = tmpp * mask;
        tmp2Mins(iq1, ib1, ib2, iac3) = tmpm * mask;
      });
  Kokkos::realloc(tmp1Plus, 0, 0, 0, 0);
  Kokkos::realloc(tmp1Mins, 0, 0, 0, 0);

  Kokkos::View<Kokkos::complex<double> ****> vPlus("vp", nq1, maxnb1, nb2,
                                                   maxnb3Plus),
      vMins("vm", nq1, maxnb1, nb2, maxnb3Mins);
  Kokkos::parallel_for(
      "vploop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, nb2, maxnb3Plus}),
      KOKKOS_LAMBDA(int iq1, int ib1, int ib2, int ib3) {
        int mask = ib1 < nb1s(iq1) && ib3 < nb3Pluss(iq1);
        Kokkos::complex<double> tmpp = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpp += tmp2Plus(iq1, ib1, ib2, iac3) * Kokkos::conj(ev3Pluss(iq1, ib3, iac3));
        }
        vPlus(iq1, ib1, ib2, ib3) = tmpp * mask;
      });

  Kokkos::parallel_for(
      "vmloop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, nb2, maxnb3Mins}),
      KOKKOS_LAMBDA(int iq1, int ib1, int ib2, int ib3) {
        int mask = ib1 < nb1s(iq1) && ib3 < nb3Minss(iq1);

        Kokkos::complex<double> tmpp = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpp += tmp2Mins(iq1, ib1, ib2, iac3) * Kokkos::conj(ev3Minss(iq1, ib3, iac3));
        }
        vMins(iq1, ib1, ib2, ib3) = tmpp * mask;
      });
  Kokkos::realloc(tmp2Plus, 0, 0, 0, 0);
  Kokkos::realloc(tmp2Mins, 0, 0, 0, 0);

  Kokkos::View<double ****> couplingPlus("cp", nq1, maxnb1, nb2, maxnb3Plus),
      couplingMins("cp", nq1, maxnb1, nb2, maxnb3Mins);
  Kokkos::parallel_for(
      "cploop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, nb2, maxnb3Plus}),
      KOKKOS_LAMBDA(int iq1, int ib1, int ib2, int ib3) {
        auto tmp = vPlus(iq1, ib1, ib2, ib3);
        couplingPlus(iq1, ib1, ib2, ib3) =
            tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
      });

  Kokkos::parallel_for(
      "cmloop",
      Kokkos::MDRangePolicy<Kokkos::Rank<4>>({0, 0, 0, 0},
                                             {nq1, maxnb1, nb2, maxnb3Mins}),
      KOKKOS_LAMBDA(int iq1, int ib1, int ib2, int ib3) {
        auto tmp = vMins(iq1, ib1, ib2, ib3);
        couplingMins(iq1, ib1, ib2, ib3) =
            tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
      });
  Kokkos::realloc(vPlus, 0, 0, 0, 0);
  Kokkos::realloc(vMins, 0, 0, 0, 0);

  // Copy result to vector of Eigen tensors
  std::vector<Eigen::Tensor<double, 3>> couplingPlus_e(nq1),
      couplingMins_e(nq1);
  auto couplingPlus_h = Kokkos::create_mirror_view(couplingPlus),
       couplingMins_h = Kokkos::create_mirror_view(couplingMins);
  Kokkos::deep_copy(couplingPlus_h, couplingPlus);
  Kokkos::deep_copy(couplingMins_h, couplingMins);
  for (int iq1 = 0; iq1 < nq1; iq1++) {
    int nb1 = nb1s_e[iq1];
    int nb3Plus = nb3Pluss_e[iq1];
    int nb3Mins = nb3Minss_e[iq1];
    couplingPlus_e[iq1] = Eigen::Tensor<double, 3>(nb1, nb2, nb3Plus);
    couplingMins_e[iq1] = Eigen::Tensor<double, 3>(nb1, nb2, nb3Mins);
    for (int ib1 = 0; ib1 < nb1; ib1++) {
      for (int ib2 = 0; ib2 < nb2; ib2++) {
        for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
          couplingPlus_e[iq1](ib1, ib2, ib3) =
              couplingPlus_h(iq1, ib1, ib2, ib3);
        }
        for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
          couplingMins_e[iq1](ib1, ib2, ib3) =
              couplingMins_h(iq1, ib1, ib2, ib3);
        }
      }
    }
  }
  return {couplingPlus_e, couplingMins_e};
}

int Interaction3Ph::estimateNumBatches(const int &nq1, const int &nb2) {
  int maxnb1 = numBands;
  int maxnb3Plus = numBands;
  int maxnb3Mins = numBands;

  // available memory is MAXMEM minus size of D3, D3cache and ev2
  double availmem = kokkosDeviceMemory->getAvailableMemory();

  // memory used by different tensors
  // Note: 16 (2*8) is the size of double (complex<double>) in bytes
  double evs = 16 * numBands * (maxnb1 + maxnb3Plus + maxnb3Mins);
  double phase = 16 * nr3;
  double tmp = 2 * 16 * numBands * numBands * numBands;
  double tmp1 = 2 * 16 * maxnb1 * numBands * numBands;
  double tmp2 = 2 * 16 * maxnb1 * nb2 * numBands;
  double v = 16 * maxnb1 * nb2 * (maxnb3Plus + maxnb3Mins);
  double c = 16 * maxnb1 * nb2 * (maxnb3Plus + maxnb3Mins);
  double maxusage =
      nq1 * (evs + std::max({phase + tmp, tmp + tmp1, tmp1 + tmp2, tmp2 + v, v + c}));

  // the number of batches needed
  int numBatches = std::ceil(maxusage / availmem);
  double totalMemory = kokkosDeviceMemory->getTotalMemory();

  if (availmem < maxusage / nq1) {
    // not enough memory to do even a single q1
    std::cerr << "total Memory = " << totalMemory / 1e9
              << "(Gb), availmem = " << availmem / 1e9
              << "(Gb), maxusage = " << maxusage / 1e9
              << "(Gb), numBatches = " << numBatches << "\n";
    Error("Insufficient memory!");
  }
  return numBatches;
}

double Interaction3Ph::getDeviceMemoryUsage() {
  double occupiedMemory = 16 * (D3_k.size() + D3PlusCached_k.size()
                                + D3MinsCached_k.size())
      + 8 * (cellPositions2_k.size() + cellPositions2_k.size()
             + weights2_k.size() + weights3_k.size());
  return occupiedMemory;
}
