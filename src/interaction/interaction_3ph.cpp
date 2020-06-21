#include "interaction_3ph.h"

long findIndexRow(Eigen::MatrixXd &cellPositions2, Eigen::Vector3d &position2) {
  long ir2 = -1;
  for (int i = 0; i < cellPositions2.cols(); i++) {
    if ((position2 - cellPositions2.col(i)).norm() == 0.) {
      ir2 = i;
      return ir2;
    }
  }
  if (ir2 == -1) {
    Error e("index not found");
  }
  return ir2;
}

// default constructor
Interaction3Ph::Interaction3Ph(Crystal &crystal_, long &numTriplets_,
                               Eigen::Tensor<double, 4> &ifc3Tensor_,
                               Eigen::Tensor<double, 3> &cellPositions_,
                               Eigen::Tensor<long, 2> &displacedAtoms_)
    : crystal(crystal_), numTriplets(numTriplets_), ifc3Tensor(ifc3Tensor_),
      cellPositions(cellPositions_), displacedAtoms(displacedAtoms_), dts(10) {

  numAtoms = crystal.getNumAtoms();
  numBands = numAtoms * 3;

  nr2 = 0;
  nr3 = 0;
  std::vector<Eigen::Vector3d> tmpCellPositions2, tmpCellPositions3;

  for (long it = 0; it < numTriplets; it++) {
    // load the position of the 2 atom in the current triplet
    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }
    // now check if this element is in the list.
    bool found2 = false;
    if (std::find(tmpCellPositions2.begin(), tmpCellPositions2.end(),
                  position2) != tmpCellPositions2.end()) {
      found2 = true;
    }
    bool found3 = false;
    if (std::find(tmpCellPositions3.begin(), tmpCellPositions3.end(),
                  position3) != tmpCellPositions3.end()) {
      found3 = true;
    }

    if (!found2) {
      tmpCellPositions2.push_back(position2);
      nr2++;
    }
    if (!found3) {
      tmpCellPositions3.push_back(position3);
      nr3++;
    }
  }

  cellPositions2 = Eigen::MatrixXd::Zero(3, nr2);
  cellPositions3 = Eigen::MatrixXd::Zero(3, nr3);
  for (int i = 0; i < nr2; i++) {
    cellPositions2.col(i) = tmpCellPositions2[i];
  }
  for (int i = 0; i < nr3; i++) {
    cellPositions3.col(i) = tmpCellPositions3[i];
  }

  D3 = Eigen::Tensor<double, 5>(numBands, numBands, numBands, nr2, nr3);
  D3.setZero();

  for (long it = 0; it < numTriplets; it++) { // sum over all triplets
    long ia1 = displacedAtoms(it, 0);
    long ia2 = displacedAtoms(it, 1);
    long ia3 = displacedAtoms(it, 2);

    Eigen::Vector3d position2, position3;
    for (int ic : {0, 1, 2}) {
      position2(ic) = cellPositions(it, 0, ic);
      position3(ic) = cellPositions(it, 1, ic);
    }

    long ir2 = findIndexRow(cellPositions2, position2);
    long ir3 = findIndexRow(cellPositions3, position3);

    for (int ic1 : {0, 1, 2}) {
      for (int ic2 : {0, 1, 2}) {
        for (int ic3 : {0, 1, 2}) {

          auto ind1 = compress2Indeces(ia1, ic1, numAtoms, 3);
          auto ind2 = compress2Indeces(ia2, ic2, numAtoms, 3);
          auto ind3 = compress2Indeces(ia3, ic3, numAtoms, 3);

          D3(ind1, ind2, ind3, ir2, ir3) = ifc3Tensor(ic3, ic2, ic1, it);
        }
      }
    }
  }
  ifc3Tensor.resize(0, 0, 0, 0);
  cellPositions.resize(0, 0, 0);
  displacedAtoms.resize(0, 0);

  // Copy everything to kokkos views

  Kokkos::realloc(cellPositions2_k, nr2, 3);
  Kokkos::realloc(cellPositions3_k, nr3, 3);

  auto cellPositions2_h = Kokkos::create_mirror_view(cellPositions2_k);
  auto cellPositions3_h = Kokkos::create_mirror_view(cellPositions3_k);
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < nr2; i++) {
      cellPositions2_h(i, j) = cellPositions2(j, i);
    }
    for (int i = 0; i < nr3; i++) {
      cellPositions3_h(i, j) = cellPositions3(j, i);
    }
  }
  Kokkos::deep_copy(cellPositions2_k, cellPositions2_h);
  Kokkos::deep_copy(cellPositions3_k, cellPositions3_h);

  ifc3Tensor.resize(0, 0, 0, 0);
  cellPositions.resize(0, 0, 0);
  displacedAtoms.resize(0, 0);

  cachedCoords << -111., -111., -111.;

  Kokkos::realloc(D3_k, numBands, numBands, numBands, nr3, nr2);
  Kokkos::realloc(D3PlusCached_k, numBands, numBands, numBands, nr3);
  Kokkos::realloc(D3MinsCached_k, numBands, numBands, numBands, nr3);
  auto D3_h = Kokkos::create_mirror_view(D3_k);
  for (int i1 = 0; i1 < numBands; i1++) {
    for (int i2 = 0; i2 < numBands; i2++) {
      for (int i3 = 0; i3 < numBands; i3++) {
        for (int i4 = 0; i4 < nr2; i4++) {
          for (int i5 = 0; i5 < nr3; i5++) {
            D3_h(i1, i2, i3, i5, i4) = D3(i1, i2, i3, i4, i5);
          }
        }
      }
    }
  }
  Kokkos::deep_copy(D3_k, D3_h);
}

// copy constructor
Interaction3Ph::Interaction3Ph(const Interaction3Ph &that)
    : crystal(that.crystal), numTriplets(that.numTriplets),
      ifc3Tensor(that.ifc3Tensor), cellPositions(that.cellPositions),
      displacedAtoms(that.displacedAtoms), tableAtCIndex1(that.tableAtCIndex1),
      tableAtCIndex2(that.tableAtCIndex2), tableAtCIndex3(that.tableAtCIndex3),
      useD3Caching(that.useD3Caching), cellPositions2(that.cellPositions2),
      cellPositions3(that.cellPositions3), D3(that.D3), nr2(that.nr2),
      nr3(that.nr3), numAtoms(that.numAtoms), numBands(that.numBands),
      cachedCoords(that.cachedCoords), D3PlusCached(that.D3PlusCached),
      D3MinsCached(that.D3MinsCached) {
  std::cout << "copy constructor called\n";
}

// assignment operator
Interaction3Ph &Interaction3Ph::operator=(const Interaction3Ph &that) {
  if (this != &that) {
    crystal = that.crystal;
    numTriplets = that.numTriplets;
    ifc3Tensor = that.ifc3Tensor;
    cellPositions = that.cellPositions;
    displacedAtoms = that.displacedAtoms;
    tableAtCIndex1 = that.tableAtCIndex1;
    tableAtCIndex2 = that.tableAtCIndex2;
    tableAtCIndex3 = that.tableAtCIndex3;

    useD3Caching = that.useD3Caching;
    cellPositions2 = that.cellPositions2;
    cellPositions3 = that.cellPositions3;
    D3 = that.D3;
    nr2 = that.nr2;
    nr3 = that.nr3;
    numAtoms = that.numAtoms;
    numBands = that.numBands;
    cachedCoords = that.cachedCoords;
    D3PlusCached = that.D3PlusCached;
    D3MinsCached = that.D3MinsCached;
  }
  return *this;
  std::cout << "assignment operator called\n";
}

template <typename A, typename B, typename C>
std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>>
Interaction3Ph::getCouplingSquared(A &state1, B &state2, C &state3Plus,
                                   C &state3Mins) {

  Eigen::Vector3d cell2Pos, cell3Pos;

  // Cartesian phonon wave vectors: q1,q2,q3
  auto q1 = state1.getCoords(Points::cartesianCoords);
  auto q2 = state2.getCoords(Points::cartesianCoords);

  Kokkos::View<double *> q1_k("q1", 3), q2_k("q2", 3);
  auto q1_h = Kokkos::create_mirror_view(q1_k),
       q2_h = Kokkos::create_mirror_view(q2_k);
  for (int i = 0; i < 3; i++) {
    q1_h(i) = q1(i);
    q2_h(i) = q2(i);
  }
  Kokkos::deep_copy(q1_k, q1_h);
  Kokkos::deep_copy(q2_k, q2_h);

  Eigen::MatrixXcd ev1_e, ev2_e, ev3Plus_e, ev3Mins_e;

  state1.getEigenvectors(ev1_e);
  state2.getEigenvectors(ev2_e);
  state3Plus.getEigenvectors(ev3Plus_e);
  state3Mins.getEigenvectors(ev3Mins_e);

  long nb1 = ev1_e.rows();
  long nb2 = ev2_e.rows();
  long nb3Plus = ev3Plus_e.rows();
  long nb3Mins = ev3Mins_e.rows();

  Kokkos::View<Kokkos::complex<double> **> ev1("ev1", nb1, nb1),
      ev2("ev2", nb2, nb2), ev3Plus("ev3p", nb3Plus, nb3Plus),
      ev3Mins("ev3m", nb3Mins, nb3Mins);

  auto e2k = [](auto e, auto k) {
    auto h = Kokkos::create_mirror_view(k);
    for (int i = 0; i < e.rows(); i++) {
      for (int j = 0; j < e.cols(); j++) {
        h(j, i) = e(i, j);
      }
    }
    Kokkos::deep_copy(k, h);
  };

  e2k(ev1_e, ev1);
  e2k(ev2_e, ev2);
  e2k(ev3Plus_e, ev3Plus);
  e2k(ev3Mins_e, ev3Mins);

  time_point t0, t1;

  Kokkos::complex<double> complexI(0.0, 1.0);

  // printf("cache test\n");
  int nr2 = this->nr2;
  int nr3 = this->nr3;
  int numBands = this->numBands;

  auto cellPositions2_k = this->cellPositions2_k;
  auto cellPositions3_k = this->cellPositions3_k;
  auto D3_k = this->D3_k;
  auto D3PlusCached_k = this->D3PlusCached_k;
  auto D3MinsCached_k = this->D3MinsCached_k;

  if (state2.getCoords(Points::cartesianCoords) != cachedCoords) {
    cachedCoords = state2.getCoords(Points::cartesianCoords);
    Kokkos::View<Kokkos::complex<double> **> phasePlus("pp", nr3, nr2),
        phaseMins("pm", nr3, nr2);

    t0 = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "phase1loop",
        Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0, 0}, {nr3, nr2}),
        KOKKOS_LAMBDA(int ir3, int ir2) {
          double argP = 0, argM = 0;
          Kokkos::complex<double> complexI(0.0, 1.0);
          for (int ic = 0; ic < 3; ic++) {
            argP += +q2_k(ic) *
                    (cellPositions2_k(ir2, ic) - cellPositions3_k(ir3, ic));
            argM += -q2_k(ic) *
                    (cellPositions2_k(ir2, ic) - cellPositions3_k(ir3, ic));
          }
          phasePlus(ir3, ir2) = Kokkos::exp(complexI * argP);
          phaseMins(ir3, ir2) = Kokkos::exp(complexI * argM);
        });
    t1 = std::chrono::steady_clock::now();
    dts[0] += t1 - t0;
    // std::cout << phasePlus(5, 8) << ", " << phaseMins(5, 8) << "\n";

    t0 = std::chrono::steady_clock::now();
    Kokkos::parallel_for(
        "D3cacheloop",
        Kokkos::MDRangePolicy<Kokkos::Rank<4>>(
            {0, 0, 0, 0}, {numBands, numBands, numBands, nr3}),
        KOKKOS_LAMBDA(int ind1, int ind2, int ind3, int ir3) {
          Kokkos::complex<double> tmpp = 0, tmpm = 0;
          for (int ir2 = 0; ir2 < nr2; ir2++) { // sum over all triplets

            // As a convention, the first primitive cell in the triplet is
            // restricted to the origin, so the phase for that cell is
            // unity.

            // printf("ind1 = %d, ind2 = %d, ind3 = %d, ir2 = %d, ir3 =
            // %d, "
            //       "nb = "
            //       "%d, nr2 = %d, nr3 = %d\n",
            //       ind1, ind2, ind3, ir2, ir3, numBands, nr2, nr3);

            tmpp += D3_k(ind1, ind2, ind3, ir3, ir2) * phasePlus(ir3, ir2);
            tmpm += D3_k(ind1, ind2, ind3, ir3, ir2) * phaseMins(ir3, ir2);
          }
          D3PlusCached_k(ind1, ind2, ind3, ir3) = tmpp;
          D3MinsCached_k(ind1, ind2, ind3, ir3) = tmpm;
        });
    t1 = std::chrono::steady_clock::now();
    dts[1] += t1 - t0;
  }

  //  std::cout << D3PlusCached_k(1, 1, 1, 1) << ", " << D3MinsCached_k(1,
  //  1, 1, 1)
  //            << "\n";
  Kokkos::View<Kokkos::complex<double> *> phasePlus("pp", nr3),
      phaseMins("pm", nr3);

  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "phase2loop", nr3, KOKKOS_LAMBDA(int ir3) {
        double argP = 0, argM = 0;
        for (int ic : {0, 1, 2}) {
          argP += -q1_k(ic) * cellPositions3_k(ir3, ic);
          argM += -q1_k(ic) * cellPositions3_k(ir3, ic);
        }
        phasePlus(ir3) = exp(complexI * argP);
        phaseMins(ir3) = exp(complexI * argM);
      });
  t1 = std::chrono::steady_clock::now();
  dts[2] += t1 - t0;

  Kokkos::View<Kokkos::complex<double> ***> tmpPlus("tmpp", numBands, numBands,
                                                    numBands),
      tmpMins("tmpm", numBands, numBands, numBands);
  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "tmploop",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {numBands, numBands, numBands}),
      KOKKOS_LAMBDA(int iac1, int iac2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int ir3 = 0; ir3 < nr3; ir3++) { // sum over all triplets
          tmpp += D3PlusCached_k(iac1, iac2, iac3, ir3) * phasePlus(ir3);
          tmpm += D3MinsCached_k(iac1, iac2, iac3, ir3) * phaseMins(ir3);
        }
        tmpPlus(iac1, iac2, iac3) = tmpp;
        tmpMins(iac1, iac2, iac3) = tmpm;
      });
  t1 = std::chrono::steady_clock::now();
  dts[3] += t1 - t0;

  // now we want to multiply
  // vPlus(ib1,ib2,ib3) = tmpPlus(iac1,iac2,iac3) * ev1(iac1,ib1)
  //          * std::conj(ev2(iac2,ib2)) * std::conj(ev3Plus(iac3,ib3));
  // vMins(ib1,ib2,ib3) += tmpMins(iac1,iac2,iac3) * ev1(iac1,ib1)
  //			* std::conj(ev2(iac2,ib2)) *
  // std::conj(ev3Mins(iac3,ib3));
  // a single loop over all 6 indices at once is easier to read, but slow
  // instead, we break it into 3 loops with 4 nested loops each
  // ( ~numBands^2 faster)

  Kokkos::View<Kokkos::complex<double> ***> tmp1Plus("t1p", nb1, numBands,
                                                     numBands),
      tmp1Mins("t1m", nb1, numBands, numBands);

  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "tmp1loop",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0},
                                             {nb1, numBands, numBands}),
      KOKKOS_LAMBDA(int ib1, int iac2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int iac1 = 0; iac1 < numBands; iac1++) {
          tmpp += tmpPlus(iac1, iac2, iac3) * ev1(ib1, iac1);
          tmpm += tmpMins(iac1, iac2, iac3) * ev1(ib1, iac1);
        }
        tmp1Plus(ib1, iac3, iac2) = tmpp;
        tmp1Mins(ib1, iac3, iac2) = tmpm;
      });
  t1 = std::chrono::steady_clock::now();
  dts[4] += t1 - t0;
  Kokkos::View<Kokkos::complex<double> ***> tmp2Plus("t2p", nb1, nb2, numBands),
      tmp2Mins("t2m", nb1, nb2, numBands);

  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "tmp2loop",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, numBands}),
      KOKKOS_LAMBDA(int ib1, int ib2, int iac3) {
        Kokkos::complex<double> tmpp = 0, tmpm = 0;
        for (int iac2 = 0; iac2 < numBands; iac2++) {
          tmpp += tmp1Plus(ib1, iac3, iac2) * ev2(ib2, iac2);
          tmpm += tmp1Mins(ib1, iac3, iac2) * Kokkos::conj(ev2(ib2, iac2));
        }
        tmp2Plus(ib1, ib2, iac3) = tmpp;
        tmp2Mins(ib1, ib2, iac3) = tmpm;
      });
  t1 = std::chrono::steady_clock::now();
  dts[5] += t1 - t0;

  Kokkos::View<Kokkos::complex<double> ***> vPlus("vp", nb1, nb2, nb3Plus),
      vMins("vm", nb1, nb2, nb3Mins);

  // the last loop is split in two because nb3 is not the same for + and -
  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "vploop",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, nb3Plus}),
      KOKKOS_LAMBDA(int ib1, int ib2, int ib3) {
        Kokkos::complex<double> tmpp = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpp += tmp2Plus(ib1, ib2, iac3) * Kokkos::conj(ev3Plus(ib3, iac3));
        }
        vPlus(ib1, ib2, ib3) = tmpp;
      });
  t1 = std::chrono::steady_clock::now();
  dts[6] += t1 - t0;

  t0 = std::chrono::steady_clock::now();
  Kokkos::parallel_for(
      "vminloop",
      Kokkos::MDRangePolicy<Kokkos::Rank<3>>({0, 0, 0}, {nb1, nb2, nb3Mins}),
      KOKKOS_LAMBDA(int ib1, int ib2, int ib3) {
        Kokkos::complex<double> tmpm = 0;
        for (int iac3 = 0; iac3 < numBands; iac3++) {
          tmpm += tmp2Mins(ib1, ib2, iac3) * Kokkos::conj(ev3Mins(ib3, iac3));
        }
        vMins(ib1, ib2, ib3) = tmpm;
      });
  t1 = std::chrono::steady_clock::now();
  dts[7] += t1 - t0;

  Eigen::Tensor<double, 3> couplingPlus(nb1, nb2, nb3Plus);
  auto vPlus_h = Kokkos::create_mirror_view(vPlus);
  Kokkos::deep_copy(vPlus_h, vPlus);
  auto norm = [](auto x) { return x.real() * x.real() + x.imag() * x.imag(); };
  // case +
  t0 = std::chrono::steady_clock::now();
  for (int ib1 = 0; ib1 < nb1; ib1++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib3 = 0; ib3 < nb3Plus; ib3++) {
        couplingPlus(ib1, ib2, ib3) = norm(vPlus_h(ib1, ib2, ib3));
      }
    }
  }
  t1 = std::chrono::steady_clock::now();
  dts[8] += t1 - t0;

  // case -
  Eigen::Tensor<double, 3> couplingMins(nb1, nb2, nb3Mins);
  auto vMins_h = Kokkos::create_mirror_view(vMins);
  Kokkos::deep_copy(vMins_h, vMins);
  for (int ib1 = 0; ib1 < nb1; ib1++) {
    for (int ib2 = 0; ib2 < nb2; ib2++) {
      for (int ib3 = 0; ib3 < nb3Mins; ib3++) {
        couplingMins(ib1, ib2, ib3) = norm(vMins_h(ib1, ib2, ib3));
      }
    }
  }
  t1 = std::chrono::steady_clock::now();
  dts[9] += t1 - t0;
  return {couplingPlus, couplingMins};
}
typedef std::tuple<Eigen::Tensor<double, 3>, Eigen::Tensor<double, 3>> eigtuple;
template eigtuple Interaction3Ph::getCouplingSquared<State, State, State>(
    State &state1, State &state2, State &state3Plus, State &state3Mins);
template eigtuple
Interaction3Ph::getCouplingSquared<State, State, DetachedState>(
    State &state1, State &state2, DetachedState &state3Plus,
    DetachedState &state3Mins);
template eigtuple
Interaction3Ph::getCouplingSquared<DetachedState, DetachedState, DetachedState>(
    DetachedState &state1, DetachedState &state2, DetachedState &state3Plus,
    DetachedState &state3Mins);

//// Function to calculate the full set of V_minus processes for a given
/// IBZ
/// mode
// void PhInteraction3Ph::calculateAllVminus(const int grid[3], const
// PhononMode &mode, 		const Eigen::MatrixXd &qFBZ,
// const Eigen::Tensor<complex<double>,3>
//&ev,const int numTriplets, 		const Eigen::Tensor<double,4>
//&ifc3Tensor, 		const
// Eigen::Tensor<double,3> &cellPositions, 		const
// Eigen::Tensor<int,2> &displacedAtoms, const CrystalInfo &crysInfo){
//
//	int iq1,iq2,iq3,ib,jb,s1,s2,s3,idim,iat;
//	int i1x,i1y,i1z,i2x,i2y,i2z,i3x,i3y,i3z;
//
//	Eigen::Tensor<complex <double>,3>
// ev1(3,crysInfo.numAtoms,crysInfo.numBranches);
// Eigen::Tensor<complex
//<double>,3> ev2(3,crysInfo.numAtoms,crysInfo.numBranches);
//	Eigen::Tensor<complex <double>,3>
// ev3(3,crysInfo.numAtoms,crysInfo.numBranches);
//
//	// Edge lengths of BZ
//	int nx = grid[0];
//	int ny = grid[1];
//	int nz = grid[2];
//	int gridSize = nx*ny*nz;
//
//	PhononTriplet interactingPhonons;
//	PhInteraction3Ph phInt;
//
//	int numAtoms = crysInfo.numAtoms;
//	int numBranches = crysInfo.numBranches;
//	Eigen::Tensor<double,4>
// Vm2(gridSize,numBranches,gridSize,numBranches);
//
//	// Grab irred phonon mode info:
//	iq1 = mode.iq; //index of wave vector in the full BZ
//	s1 = mode.s; //branch
//
//	//Demux 1st phonon wave vector
//	//!!WARNING: For testing purposes using ShengBTE ordering!!!
//	i1x = iq1%nx;
//	i1y = (iq1/nx)%ny;
//	i1z = iq1/nx/ny;
//
//	interactingPhonons.s1 = s1;
//	interactingPhonons.iq1 = iq1;
//	for(idim = 0; idim < 3; idim++){
//		for(iat = 0; iat < numAtoms; iat++){
//			for(ib = 0; ib < numBranches; ib++){
//				ev1(idim,iat,ib) =
// ev(iq1,idim+3*iat,ib);
//			}
//		}
//	}
//	interactingPhonons.ev1 = ev1;
//
//	int count = 0;
//
//	// Loop over all 2nd phonon wave vectors
//	for(i2z = 0; i2z < nz; i2z++){
//		for(i2y = 0; i2y < ny; i2y++){
//			for(i2x = 0; i2x < nx; i2x++){
//				//Muxed index of 2nd phonon wave vector
//				//!!WARNING: For testing purposes using
// ShengBTE
// ordering!!! 				iq2 = (i2z*ny + i2y)*nx + i2x;
//
//				interactingPhonons.iq2 = iq2;
//
//				for(idim = 0; idim < 3; idim++){
//					for(iat = 0; iat < numAtoms;
// iat++){ 						for(ib = 0; ib <
// numBranches;
// ib++){
// ev2(idim,iat,ib) = ev(iq2,idim+3*iat,ib);
//						}
//					}
//				}
//				interactingPhonons.ev2 = ev2;
//
//				// Third phonon wave vector (Umklapped,
// if
// needed) 				i3x = (i1x - i2x + nx)%nx;
// i3y = (i1y - i2y + ny)%ny; 				i3z = (i1z - i2z
//+ nz)%nz;
//				//!!WARNING: For testing purposes using
// ShengBTE
// ordering!!! 				iq3 = (i3z*ny + i3y)*nx + i3x;
//
//				interactingPhonons.iq3 = iq3;
//				for(idim = 0; idim < 3; idim++){
//					for(iat = 0; iat < numAtoms;
// iat++){ 						for(ib = 0; ib <
// numBranches;
// ib++){
// ev3(idim,iat,ib) = ev(iq3,idim+3*iat,ib);
//						}
//					}
//				}
//				interactingPhonons.ev3 = ev3;
//
//				// Sum over 2nd phonon branches
//				for(ib = 0; ib < numBranches; ib++){
//					interactingPhonons.s2 = ib;
//					// Sum over 3rd phonon branches
//					for(jb = 0; jb < numBranches;
// jb++){ 						interactingPhonons.s3 =
// jb;
//
//						// Call calculateSingleV
//						Vm2(iq2,ib,iq3,jb) =
// phInt.calculateSingleV(interactingPhonons, qFBZ, numTriplets,
// ifc3Tensor,
//								cellPositions,
// displacedAtoms, crysInfo, '-');
//
//						//cout << iq2+1 << " "
//<< ib+1
//<<
//"
//"
//<< iq3+1 << " " << jb+1 << " " << Vm2(iq2,ib,iq3,jb) << "\n";
//					}
//				}
//			}
//		}
//		cout << ++count*2304.0/18432.0*100 << " % done.\n";
//	}
//
//
//	// Write to disk
//	string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
//	ofstream outFile;
//	//outFile.open(fileName, ios::out | ios::trunc | ios::binary);
//	outFile.open(fileName,ios::trunc);
//	for(iq2 = 0; iq2 < gridSize; iq2++){
//		for(iq3 = 0; iq3 < gridSize; iq3++){
//			for(ib = 0; ib < numBranches; ib++){
//				for(jb = 0; jb < numBranches; jb++){
//					outFile << Vm2(iq2,ib,iq3,jb) <<
//"\n";
//				}
//			}
//		}
//	}
//	outFile.close();
//}
//
////Transition probabilities for a given irreducible phonon mode
// void PhInteraction3Ph::calculateAllW(const double T,const int grid[3],
// const PhononMode &mode, 		const Eigen::MatrixXi
// &indexMesh,const CrystalInfo
//&crysInfo, 		const Eigen::MatrixXd omega,const TetraData
// tetra){
//
//	int
// iq1,iq2,iq3,iq3Plus,iq3Minus,s1,s2,s3,iDim,plusProcessCount,minusProcessCount;
//	int numBranches = crysInfo.numBranches;
//	int nq = grid[0]*grid[1]*grid[2];
//
//	double
// omega1,omega2,omega3Plus,omega3Minus,n01,n02,n03Plus,n03Minus,Vplus2,Vminus2,
//	Wplus,Wminus,tetWeightPlus,tetWeightMinus;
//	const double a = M_PI*hbar/4.0*5.60626442*1.0e30;
//
//	Eigen::Vector3i q1,q2,q3Plus,q3Minus;
//
//	// Grab irred phonon mode info:
//	iq1 = mode.iq; //index of wave vector in the full BZ
//	s1 = mode.s; //branch
//	omega1 = omega(iq1,s1); //irred mode energy
//	q1 = indexMesh.row(iq1);
//
//	//Read full set of V-(q2,s2;q3,s3) for given mode from file
//	Eigen::Tensor<double,4> Vm2(nq,numBranches,nq,numBranches);
//	string fileName = "Vm2.iq"+to_string(iq1)+".s"+to_string(s1);
//	ifstream inFile;
//	inFile.open(fileName);
//	for(iq2 = 0; iq2 < nq; iq2++){
//		for(iq3 = 0; iq3 < nq; iq3++){
//			for(s2 = 0; s2 < numBranches; s2++){
//				for(s3 = 0; s3 < numBranches; s3++){
//					inFile >> Vm2(iq2,s2,iq3,s3);
//				}
//			}
//		}
//	}
//	inFile.close();
//
//	//Open W+, W-, and process counter files
//	string WplusFileName =
//"Wp2.iq"+to_string(iq1)+".s"+to_string(s1); 	string WminusFileName =
//"Wm2.iq"+to_string(iq1)+".s"+to_string(s1); 	ofstream WplusFile,
// WminusFile; 	WplusFile.open(WplusFileName, ios::trunc);
//	WminusFile.open(WminusFileName, ios::trunc);
//
//	plusProcessCount = 0;
//	minusProcessCount = 0;
//	if(omega1 > 0){ //skip zero energy phonon
//		//Bose distribution for first phonon mode
//		n01 = bose(omega1,T);
//
//		//Sum over second phonon wave vector in full BZ
//		for(iq2 = 0; iq2 < nq; iq2++){
//			q2 = indexMesh.row(iq2);
//
//			//Get final phonon wave vector location
//			//modulo reciprocal lattice vector
//			for(iDim = 0; iDim < 3; iDim++){
//				//plus process
//				q3Plus(iDim) =
//(q1(iDim)+q2(iDim))%grid[iDim];
//				//minus process
//				q3Minus(iDim) =
//(q1(iDim)-q2(iDim)+grid[iDim])%grid[iDim];
//			}
//			//!!WARNING: For testing purposes using ShengBTE
// ordering!!! 			iq3Plus = (q3Plus(2)*grid[1] +
// q3Plus(1))*grid[0]
// + q3Plus(0); 			iq3Minus = (q3Minus(2)*grid[1] +
// q3Minus(1))*grid[0] + q3Minus(0);
//
//			//Sum over second phonon branches
//			for(s2 = 0; s2 < numBranches; s2++){
//				omega2 = omega(iq2,s2); //second phonon
// ang.
// freq.
//
//				if(omega2 > 0){ //skip zero energy
// phonon
//
//					//Bose distribution for second
// phonon
// mode 					n02 = bose(omega2,T);
//
//					//Sum over final phonon branches
//					for(s3 = 0; s3 < numBranches;
// s3++){ 						omega3Plus =
// omega(iq3Plus,s3);
////third
// phonon(+) ang. freq. omega3Minus = omega(iq3Minus,s3); //third
// phonon(-) ang. freq.
//
//						//Bose distribution for
// final
// phonon mode 						n03Plus =
// bose(omega3Plus,T); 						n03Minus
// = bose(omega3Minus,T);
//
//						//Calculate tetrahedron
// weight
// for plus
// and minus processes tetWeightPlus =
// fillTetsWeights(omega3Plus-omega1,s2,iq2,tetra);
// tetWeightMinus = fillTetsWeights(omega1-omega3Minus,s2,iq2,tetra);
//
//						//Plus processes
//						if(tetWeightPlus > 0 &&
// omega3Plus > 0){
//							//Increase
// processes
// counter plusProcessCount++;
//
//							//Time reverse
// second
// phonon to get Vplus2 Vplus2
// =
// Vm2(timeReverse(iq2,grid),s2,iq3Plus,s3);
//
//							//Calculatate
// transition
// probability W+ 							Wplus
// = a*(n02-n03Plus)*Vplus2*tetWeightPlus
//									/(omega1*omega2*omega3Plus);
////THz
//
//							//Write plus
// process
// info to
// file 							WplusFile << iq2
// <<
// "
// "
// << s2
// <<
// "
// "
// << iq3Plus << " " << s3 << " "
//									<<
// Wplus
//<<
//"\n";
//						}
//
//						//Minus processes
//						if(tetWeightMinus > 0 &&
// omega3Minus > 0){
//							//Increase
// processes
// counter minusProcessCount++;
//
//							Vminus2 =
// Vm2(iq2,s2,iq3Minus,s3);
//
//							//Calculatate
// transition
// probability W- 							Wminus
// = a*(n02+n03Minus+1.0)*Vminus2*tetWeightMinus
//									/(omega1*omega2*omega3Minus);
////THz
//
//							//Write minus
// process
// info to
// disk 							WminusFile <<
// iq2
// <<
// "
// "
// << s2
// <<
// "
// "
// << iq3Minus << " " << s3 << " "
//									<<
// Wminus
//<<
//"\n";
//						}
//					}//s3
//				}//zero of second phonon
//			}//s2
//		}//iq2
//	}//zero of first phonon
//	WplusFile.close();
//	WminusFile.close();
//
//	//Write total number of plus and minus processes to disk
//	string counterFileName =
//"WCounter.iq"+to_string(iq1)+".s"+to_string(s1); 	ofstream
// counterFile; 	counterFile.open(counterFileName, ios::trunc);
// counterFile
//<< plusProcessCount << "\n"; 	counterFile << minusProcessCount << "\n";
//	counterFile.close();
//}
