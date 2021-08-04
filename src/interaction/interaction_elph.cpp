#include "interaction_elph.h"
#include <Kokkos_Core.hpp>
#include <fstream>

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

// default constructor
InteractionElPhWan::InteractionElPhWan(
    Crystal &crystal_,
    const Eigen::Tensor<std::complex<double>, 5> &couplingWannier_,
    const Eigen::MatrixXd &elBravaisVectors_,
    const Eigen::VectorXd &elBravaisVectorsDegeneracies_,
    const Eigen::MatrixXd &phBravaisVectors_,
    const Eigen::VectorXd &phBravaisVectorsDegeneracies_, PhononH0 *phononH0_)
    : crystal(crystal_), phononH0(phononH0_) {

  couplingWannier = couplingWannier_;
  elBravaisVectors = elBravaisVectors_;
  elBravaisVectorsDegeneracies = elBravaisVectorsDegeneracies_;
  phBravaisVectors = phBravaisVectors_;
  phBravaisVectorsDegeneracies = phBravaisVectorsDegeneracies_;

  numElBands = couplingWannier.dimension(0);
  numPhBands = couplingWannier.dimension(2);
  numPhBravaisVectors = couplingWannier.dimension(3);
  numElBravaisVectors = couplingWannier.dimension(4);
  cachedK1.setZero();

  usePolarCorrection = false;
  if (phononH0 != nullptr) {
    Eigen::Matrix3d epsilon = phononH0->getDielectricMatrix();
    if (epsilon.squaredNorm() > 1.0e-10) { // i.e. if epsilon wasn't computed
      if (crystal.getNumSpecies() > 1) {   // otherwise polar correction = 0
        usePolarCorrection = true;
      }
    }
  }
}

InteractionElPhWan::InteractionElPhWan(Crystal &crystal_) : crystal(crystal_) {}

// copy constructor
InteractionElPhWan::InteractionElPhWan(const InteractionElPhWan &that)
    : crystal(that.crystal), phononH0(that.phononH0),
      couplingWannier(that.couplingWannier),
      couplingWannier_k(that.couplingWannier_k),
      elBravaisVectors(that.elBravaisVectors),
      elBravaisVectorsDegeneracies(that.elBravaisVectorsDegeneracies),
      phBravaisVectors(that.phBravaisVectors),
      phBravaisVectorsDegeneracies(that.phBravaisVectorsDegeneracies),
      numPhBands(that.numPhBands), numElBands(that.numElBands),
      numElBravaisVectors(that.numElBravaisVectors),
      numPhBravaisVectors(that.numPhBravaisVectors),
      cacheCoupling(that.cacheCoupling), elPhCached(that.elPhCached),
      cachedK1(that.cachedK1), usePolarCorrection(that.usePolarCorrection) {}

// assignment operator
InteractionElPhWan &
InteractionElPhWan::operator=(const InteractionElPhWan &that) {
  if (this != &that) {
    crystal = that.crystal;
    phononH0 = that.phononH0;
    couplingWannier = that.couplingWannier;
    couplingWannier_k = that.couplingWannier_k;
    elBravaisVectors = that.elBravaisVectors;
    elBravaisVectorsDegeneracies = that.elBravaisVectorsDegeneracies;
    phBravaisVectors = that.phBravaisVectors;
    phBravaisVectorsDegeneracies = that.phBravaisVectorsDegeneracies;
    numPhBands = that.numPhBands;
    numElBands = that.numElBands;
    numElBravaisVectors = that.numElBravaisVectors;
    numPhBravaisVectors = that.numPhBravaisVectors;
    cacheCoupling = that.cacheCoupling;
    elPhCached = that.elPhCached;
    cachedK1 = that.cachedK1;
    usePolarCorrection = that.usePolarCorrection;
  }
  return *this;
}

Eigen::Tensor<double, 3>
InteractionElPhWan::getCouplingSquared(const int &ik2) {
  return cacheCoupling[ik2];
}

Eigen::Tensor<std::complex<double>, 3> InteractionElPhWan::getPolarCorrection(
    const Eigen::Vector3d &q3, const Eigen::MatrixXcd &ev1,
    const Eigen::MatrixXcd &ev2, const Eigen::MatrixXcd &ev3) {
  // doi:10.1103/physrevlett.115.176401, Eq. 4, is implemented here

  // gather variables
  double volume = crystal.getVolumeUnitCell();
  Eigen::Matrix3d reciprocalUnitCell = crystal.getReciprocalUnitCell();
  Eigen::Matrix3d epsilon = phononH0->getDielectricMatrix();
  Eigen::Tensor<double, 3> bornCharges = phononH0->getBornCharges();
  // must be in Bohr
  Eigen::MatrixXd atomicPositions = crystal.getAtomicPositions();
  Eigen::Vector3i qCoarseMesh = phononH0->getCoarseGrid();

  return getPolarCorrectionStatic(q3, ev1, ev2, ev3, volume, reciprocalUnitCell,
                                  epsilon, bornCharges, atomicPositions,
                                  qCoarseMesh);
}

Eigen::Tensor<std::complex<double>, 3>
InteractionElPhWan::getPolarCorrectionStatic(
    const Eigen::Vector3d &q3, const Eigen::MatrixXcd &ev1,
    const Eigen::MatrixXcd &ev2, const Eigen::MatrixXcd &ev3,
    const double &volume, const Eigen::Matrix3d &reciprocalUnitCell,
    const Eigen::Matrix3d &epsilon, const Eigen::Tensor<double, 3> &bornCharges,
    const Eigen::MatrixXd &atomicPositions,
    const Eigen::Vector3i &qCoarseMesh) {
  // doi:10.1103/physrevlett.115.176401, Eq. 4, is implemented here

  int numAtoms = atomicPositions.rows();

  // overlap = <U^+_{b2 k+q}|U_{b1 k}>
  //         = <psi_{b2 k+q}|e^{i(q+G)r}|psi_{b1 k}>
  Eigen::MatrixXcd overlap = ev2.adjoint() * ev1; // matrix size (nb2,nb1)
  overlap = overlap.transpose();                  // matrix size (nb1,nb2)

  // auxiliary terms
  double gMax = 14.;
  double e2 = 2.; // = e^2/4/Pi/eps_0 in atomic units
  std::complex<double> factor = e2 * fourPi / volume * complexI;

  // build a list of (q+G) vectors
  std::vector<Eigen::Vector3d> gVectors; // here we insert all (q+G)
  for (int m1 = -qCoarseMesh(0); m1 <= qCoarseMesh(0); m1++) {
    for (int m2 = -qCoarseMesh(1); m2 <= qCoarseMesh(1); m2++) {
      for (int m3 = -qCoarseMesh(2); m3 <= qCoarseMesh(2); m3++) {
        Eigen::Vector3d gVector;
        gVector << m1, m2, m3;
        gVector = reciprocalUnitCell * gVector;
        gVector += q3;
        gVectors.push_back(gVector);
      }
    }
  }

  int numPhBands = ev3.rows();
  Eigen::VectorXcd x(numPhBands);
  x.setZero();
  for (Eigen::Vector3d gVector : gVectors) {
    double qEq = gVector.transpose() * epsilon * gVector;
    if (qEq > 0. && qEq / 4. < gMax) {
      std::complex<double> factor2 = factor * exp(-qEq / 4.) / qEq;
      for (int iAt = 0; iAt < numAtoms; iAt++) {
        double arg = -gVector.dot(atomicPositions.row(iAt));
        std::complex<double> phase = {cos(arg), sin(arg)};
        std::complex<double> factor3 = factor2 * phase;
        for (int iPol : {0, 1, 2}) {
          double gqDotZ = gVector(0) * bornCharges(iAt, 0, iPol) +
                          gVector(1) * bornCharges(iAt, 1, iPol) +
                          gVector(2) * bornCharges(iAt, 2, iPol);
          int k = PhononH0::getIndexEigenvector(iAt, iPol, numAtoms);
          for (int ib3 = 0; ib3 < numPhBands; ib3++) {
            x(ib3) += factor3 * gqDotZ * ev3(k, ib3);
          }
        }
      }
    }
  }

  Eigen::Tensor<std::complex<double>, 3> v(overlap.rows(), overlap.cols(),
                                           numPhBands);
  v.setZero();
  for (int ib3 = 0; ib3 < numPhBands; ib3++) {
    for (int i = 0; i < overlap.rows(); i++) {
      for (int j = 0; j < overlap.cols(); j++) {
        v(i, j, ib3) += x(ib3) * overlap(i, j);
      }
    }
  }
  return v;
}

// Forward declare these helper functions, as it reads nicely to have
// the general parse function first
InteractionElPhWan parseHDF5(Context &context, Crystal &crystal,
                             PhononH0 *phononH0_);
InteractionElPhWan parseNoHDF5(Context &context, Crystal &crystal,
                               PhononH0 *phononH0_);

// General parse function
InteractionElPhWan InteractionElPhWan::parse(Context &context, Crystal &crystal,
                                             PhononH0 *phononH0_) {
  if (mpi->mpiHead()) {
    std::cout << "\n";
    std::cout << "Started parsing of el-ph interaction." << std::endl;
  }
#ifdef HDF5_AVAIL
  auto output = parseHDF5(context, crystal, phononH0_);
#else
  auto output = parseNoHDF5(context, crystal, phononH0_);
#endif

  if (mpi->mpiHead()) {
    std::cout << "Finished parsing of el-ph interaction." << std::endl;
  }

  return output;
}
// specific parse function for the case where there is no
// HDF5 available
InteractionElPhWan parseNoHDF5(Context &context, Crystal &crystal,
                               PhononH0 *phononH0_) {

  std::string fileName = context.getElphFileName();

  int numElectrons, numSpin;
  int numElBands, numElBravaisVectors, numPhBands, numPhBravaisVectors;
  Eigen::MatrixXd phBravaisVectors_, elBravaisVectors_;
  Eigen::VectorXd phBravaisVectorsDegeneracies_, elBravaisVectorsDegeneracies_;
  Eigen::Tensor<std::complex<double>, 5> couplingWannier_;

  // Open ElPh file
  if (mpi->mpiHeadPool()) {
    std::ifstream infile(fileName);
    if (not infile.is_open()) {
      Error("ElPh file not found");
    }

    // Read the bravais lattice vectors info for q mesh.
    infile >> numElectrons >> numSpin;

    int kx, ky, kz;
    int qx, qy, qz;
    infile >> kx >> ky >> kz;
    infile >> qx >> qy >> qz;

    int iCart;

    infile >> iCart >> numPhBravaisVectors;
    phBravaisVectors_.resize(3, numPhBravaisVectors);
    phBravaisVectorsDegeneracies_.resize(numPhBravaisVectors);
    phBravaisVectors_.setZero();
    phBravaisVectorsDegeneracies_.setZero();
    for (int i : {0, 1, 2}) {
      for (int j = 0; j < numPhBravaisVectors; j++) {
        infile >> phBravaisVectors_(i, j);
      }
    }
    for (int i = 0; i < numPhBravaisVectors; i++) {
      infile >> phBravaisVectorsDegeneracies_(i);
    }

    int totalNumElBravaisVectors;
    infile >> iCart >> totalNumElBravaisVectors;

    auto localElVectors = mpi->divideWorkIter(totalNumElBravaisVectors, mpi->intraPoolComm);
    numElBravaisVectors = int(localElVectors.size());

    elBravaisVectors_.resize(3, numElBravaisVectors);
    elBravaisVectorsDegeneracies_.resize(numElBravaisVectors);
    elBravaisVectors_.setZero();
    elBravaisVectorsDegeneracies_.setZero();
    for (int i : {0, 1, 2}) {
      for (int j = 0; j < totalNumElBravaisVectors; j++) {
        double x;
        infile >> x;
        if (std::find(localElVectors.begin(),localElVectors.end(), j) != localElVectors.end() ) {
          int localIrE = j - localElVectors[0];
          elBravaisVectors_(i, localIrE) = x;
        }
      }
    }
    for (int i = 0; i < totalNumElBravaisVectors; i++) {
      double x;
      infile >> x;
      if (std::find(localElVectors.begin(),localElVectors.end(), i) != localElVectors.end() ) {
        int localIrE = i - localElVectors[0];
        elBravaisVectorsDegeneracies_(localIrE) = x;
      }
    }
    std::string line;
    std::getline(infile, line);

    // Read real space matrix elements for el-ph coupling
    int tmpI;
    infile >> numElBands >> tmpI >> numPhBands >> tmpI >> tmpI;

    // user info about memory
    {
      std::complex<double> cx;
      double x = numElBands * numElBands * numPhBands * numPhBravaisVectors *
                 numElBravaisVectors / pow(1024., 3) * sizeof(cx);
      std::cout << "Allocating " << x
                << " (GB) (per MPI process) for the el-ph coupling matrix."
                << std::endl;
    }

    couplingWannier_.resize(numElBands, numElBands, numPhBands,
                            numPhBravaisVectors, numElBravaisVectors);
    couplingWannier_.setZero();
    double re, im;
    for (int i5 = 0; i5 < totalNumElBravaisVectors; i5++) {
      int localIrE = -1;
      if (std::find(localElVectors.begin(),localElVectors.end(), i5) != localElVectors.end() ) {
        localIrE = i5 - localElVectors[0];
      }
      for (int i4 = 0; i4 < numPhBravaisVectors; i4++) {
        for (int i3 = 0; i3 < numPhBands; i3++) {
          for (int i2 = 0; i2 < numElBands; i2++) {
            for (int i1 = 0; i1 < numElBands; i1++) {
              infile >> re >> im;
              // note: in qe2Phoebe, the first index is on k+q bands,
              // and the second is on the bands of k. Here I invert them
              // similarly, in qe2Phoebe I inverted the order of R_el and R_ph
              if (localIrE >= 0) {
                couplingWannier_(i1, i2, i3, i4, localIrE) = {re, im};
              }
            }
          }
        }
      }
    }
  } // mpiHead done reading file

  mpi->bcast(&numElectrons);
  mpi->bcast(&numSpin);

  mpi->bcast(&numElBands);
  mpi->bcast(&numPhBands);
  mpi->bcast(&numElBravaisVectors, mpi->interPoolComm);
  mpi->bcast(&numPhBravaisVectors);

  if (numSpin == 2) {
    Error("Spin is not currently supported");
  }
  context.setNumOccupiedStates(numElectrons);

  if (!mpi->mpiHeadPool()) { // head already allocated these
    phBravaisVectors_.resize(3, numPhBravaisVectors);
    phBravaisVectorsDegeneracies_.resize(numPhBravaisVectors);
    elBravaisVectors_.resize(3, numElBravaisVectors);
    elBravaisVectorsDegeneracies_.resize(numElBravaisVectors);
    couplingWannier_.resize(numElBands, numElBands, numPhBands,
                            numPhBravaisVectors, numElBravaisVectors);
    phBravaisVectors_.setZero();
    phBravaisVectorsDegeneracies_.setZero();
    elBravaisVectors_.setZero();
    elBravaisVectorsDegeneracies_.setZero();
    couplingWannier_.setZero();
  }
  mpi->bcast(&elBravaisVectors_, mpi->interPoolComm);
  mpi->bcast(&elBravaisVectorsDegeneracies_, mpi->interPoolComm);
  mpi->bcast(&phBravaisVectors_);
  mpi->bcast(&phBravaisVectorsDegeneracies_);
  mpi->bcast(&couplingWannier_, mpi->interPoolComm);

  InteractionElPhWan output(crystal, couplingWannier_, elBravaisVectors_,
                            elBravaisVectorsDegeneracies_, phBravaisVectors_,
                            phBravaisVectorsDegeneracies_, phononH0_);
  return output;
}

#ifdef HDF5_AVAIL
// specific parse function for the case where parallel HDF5 is available
InteractionElPhWan parseHDF5(Context &context, Crystal &crystal,
                             PhononH0 *phononH0_) {

  std::string fileName = context.getElphFileName();

  int numElectrons, numSpin;
  int numElBands, numElBravaisVectors, numPhBands, numPhBravaisVectors;
  Eigen::MatrixXd phBravaisVectors_, elBravaisVectors_;
  Eigen::VectorXd phBravaisVectorsDegeneracies_, elBravaisVectorsDegeneracies_;
  Eigen::Tensor<std::complex<double>, 5> couplingWannier_;

  // check for existence of file
  {
    std::ifstream infile(fileName);
    if (not infile.is_open()) {
      Error("Required electron-phonon file ***.phoebe.elph.hdf5 "
            "not found at " +
            fileName + " .");
    }
  }

  try {
    // Use MPI head only to read in the small data structures
    // then distribute them below this
    if (mpi->mpiHead()) {
      // need to open the files differently if MPI is available or not
      // NOTE: do not remove the braces inside this if -- the file must
      // go out of scope, so that it can be reopened for parallel
      // read in the next block.
      {
        // Open the HDF5 ElPh file
        HighFive::File file(fileName, HighFive::File::ReadOnly);

        // read in the number of electrons and the spin
        HighFive::DataSet dnelec = file.getDataSet("/numElectrons");
        HighFive::DataSet dnspin = file.getDataSet("/numSpin");
        dnelec.read(numElectrons);
        dnspin.read(numSpin);

        // read in the number of phonon and electron bands
        HighFive::DataSet dnElBands = file.getDataSet("/numElBands");
        HighFive::DataSet dnModes = file.getDataSet("/numPhModes");
        dnElBands.read(numElBands);
        dnModes.read(numPhBands);

        // read in bravais lattice vectors
        HighFive::DataSet dphbravais = file.getDataSet("/phBravaisVectors");
        HighFive::DataSet delbravais = file.getDataSet("/elBravaisVectors");
        dphbravais.read(phBravaisVectors_);
        delbravais.read(elBravaisVectors_);
        numElBravaisVectors = elBravaisVectors_.cols();
        numPhBravaisVectors = phBravaisVectors_.cols();

        // Read in electron and phonon degeneracies
        HighFive::DataSet dphDegeneracies = file.getDataSet("/phDegeneracies");
        HighFive::DataSet delDegeneracies = file.getDataSet("/elDegeneracies");
        dphDegeneracies.read(phBravaisVectorsDegeneracies_);
        delDegeneracies.read(elBravaisVectorsDegeneracies_);
      }
    }
    // bcast to all MPI processes
    mpi->bcast(&numElectrons);
    mpi->bcast(&numSpin);
    mpi->bcast(&numElBands);
    mpi->bcast(&numPhBands);
    mpi->bcast(&numElBravaisVectors);
    mpi->bcast(&numPhBravaisVectors);

    if (numSpin == 2) {
      Error("Spin is not currently supported");
    }
    context.setNumOccupiedStates(numElectrons);

    if (!mpi->mpiHead()) { // head already allocated these
      phBravaisVectors_.resize(3, numElBravaisVectors);
      phBravaisVectorsDegeneracies_.resize(numElBravaisVectors);
      elBravaisVectors_.resize(3, numElBravaisVectors);
      elBravaisVectorsDegeneracies_.resize(numElBravaisVectors);
      couplingWannier_.resize(numElBands, numElBands, numPhBands,
                              numPhBravaisVectors, numElBravaisVectors);
      phBravaisVectors_.setZero();
      phBravaisVectorsDegeneracies_.setZero();
      elBravaisVectors_.setZero();
      elBravaisVectorsDegeneracies_.setZero();
    }
    mpi->bcast(&elBravaisVectors_);
    mpi->bcast(&elBravaisVectorsDegeneracies_);
    mpi->bcast(&phBravaisVectors_);
    mpi->bcast(&phBravaisVectorsDegeneracies_);

    // Define the eph matrix element containers
    size_t totElems = numElBands * numElBands * numPhBands *
                      numPhBravaisVectors * numElBravaisVectors;

    // user info about memory
    {
      std::complex<double> cx;
      double x = totElems / pow(1024., 3) * sizeof(cx);
      if (mpi->mpiHead()) {
        std::cout << "Allocating " << x
                  << " (GB) (per MPI process) for the el-ph coupling matrix."
                  << std::endl;
      }
    }

    couplingWannier_.resize(numElBands, numElBands, numPhBands,
                            numPhBravaisVectors, numElBravaisVectors);
    couplingWannier_.setZero();

// Regular parallel read
#if defined(MPI_AVAIL) && !defined(HDF5_SERIAL)

    // Reopen the HDF5 ElPh file for parallel read of eph matrix elements
    HighFive::File file(
        fileName, HighFive::File::ReadOnly,
        HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

    // get the start and stop points of elements to be written by this process
    std::vector<int> workDivs = mpi->divideWork(totElems);
    size_t localElems = workDivs[1] - workDivs[0];

    // Set up buffer to be filled from hdf5
    std::vector<std::complex<double>> gWanSlice(localElems);
    // Set up buffer to receive full matrix data
    std::vector<std::complex<double>> gWanFlat(totElems);

    // Set up dataset for gWannier
    HighFive::DataSet dgWannier = file.getDataSet("/gWannier");
    // Read in the elements for this process
    dgWannier.select({0, size_t(workDivs[0])}, {1, localElems})
        .read(gWanSlice);

    // Gather the elements read in by each process
    mpi->allGatherv(&gWanSlice, &gWanFlat);

    // Map the flattened matrix back to tensor structure
    Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 5>> gWanTemp(
        gWanFlat.data(), numElBands, numElBands, numPhBands,
        numPhBravaisVectors, numElBravaisVectors);
    couplingWannier_ = gWanTemp;

#else
    // Reopen serial version, either because MPI does not exist
    // or because we forced HDF5 to run in serial.

    // Set up buffer to receive full matrix data
    std::vector<std::complex<double>> gWanFlat(totElems);

    if (mpi->mpiHead()) {
      HighFive::File file(fileName, HighFive::File::ReadOnly);

      // Set up dataset for gWannier
      HighFive::DataSet dgWannier = file.getDataSet("/gWannier");
      // Read in the elements for this process
      dgWannier.read(gWanFlat);
    }
    mpi->bcast(&gWanFlat);

    // Map the flattened matrix back to tensor structure
    Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 5>> gWanTemp(
        gWanFlat.data(), numElBands, numElBands, numPhBands,
        numPhBravaisVectors, numElBravaisVectors);
    couplingWannier_ = gWanTemp;

#endif

  } catch (std::exception &error) {
    Error("Issue reading elph Wannier represenation from hdf5.");
  }

  InteractionElPhWan output(crystal, couplingWannier_, elBravaisVectors_,
                            elBravaisVectorsDegeneracies_, phBravaisVectors_,
                            phBravaisVectorsDegeneracies_, phononH0_);
  return output;
}
#endif

void InteractionElPhWan::calcCouplingSquared(
    const Eigen::MatrixXcd &eigvec1,
    const std::vector<Eigen::MatrixXcd> &eigvecs2,
    const std::vector<Eigen::MatrixXcd> &eigvecs3, const Eigen::Vector3d &k1C,
    const std::vector<Eigen::Vector3d> &k2Cs,
    const std::vector<Eigen::Vector3d> &q3Cs) {
  (void)k2Cs;
  int numWannier = numElBands;
  int nb1 = eigvec1.cols();
  int numLoops = eigvecs2.size();
  cacheCoupling.resize(0);
  cacheCoupling.resize(numLoops);
  Kokkos::complex<double> complexI(0.0, 1.0);

  auto elPhCached = this->elPhCached;
  int numPhBands = this->numPhBands;
  numElBands = this->numElBands;
  numElBravaisVectors = this->numElBravaisVectors;
  numPhBravaisVectors = this->numPhBravaisVectors;

  if (k1C != cachedK1 || elPhCached.size() == 0) {
    cachedK1 = k1C;

    for (int iPool=0; iPool<mpi->getSize(mpi->intraPoolComm); iPool++) {

      // here we decide which k1 and ev1 to use now.
      Eigen::Vector3d poolK1C = Eigen::Vector3d::Zero();
      int poolNb1 = 0;
      if (iPool == mpi->getRank(mpi->intraPoolComm)) {
        poolNb1 = nb1;
      }
      mpi->allReduceSum(&poolNb1, mpi->intraPoolComm);

      Eigen::MatrixXcd poolEigvec1 = Eigen::MatrixXcd::Zero(poolNb1, numWannier);
      if (iPool == mpi->getRank(mpi->intraPoolComm)) {
        poolK1C = k1C;
        poolEigvec1 = eigvec1;
      }
      mpi->allReduceSum(&poolK1C, mpi->intraPoolComm);
      mpi->allReduceSum(&poolEigvec1, mpi->intraPoolComm);

      ComplexView4D g1(Kokkos::ViewAllocateWithoutInitializing("g1"), numPhBravaisVectors, numPhBands, numWannier,
                       numWannier);

      ComplexView2D eigvec1_k("ev1", poolNb1, numWannier);
      {
        HostComplexView2D eigvec1_h((Kokkos::complex<double>*) poolEigvec1.data(), poolNb1, numWannier);
        Kokkos::deep_copy(eigvec1_k, eigvec1_h);
      }

      ComplexView1D phases_k("phases", numElBravaisVectors);
      auto phases_h = Kokkos::create_mirror_view(phases_k);
      for (int irE = 0; irE < numElBravaisVectors; irE++) {
        double arg = poolK1C.dot(elBravaisVectors.col(irE));
        phases_h(irE) =
            exp(complexI * arg) / double(elBravaisVectorsDegeneracies(irE));
      }
      Kokkos::deep_copy(phases_k, phases_h);

      if(couplingWannier_k.extent(0)==0){
        HostComplexView5D couplingWannier_h((Kokkos::complex<double>*) couplingWannier.data(),
            numElBravaisVectors, numPhBravaisVectors,
            numPhBands, numWannier, numWannier);
        couplingWannier_k = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), couplingWannier_h);
        Kokkos::deep_copy(couplingWannier_k, couplingWannier_h);
      }
      ComplexView5D couplingWannier_k = this->couplingWannier_k;
#ifdef KOKKOS_ENABLE_CUDA
    Kokkos::parallel_for(
        "g1",
        Range4D({0, 0, 0, 0},
                {numPhBravaisVectors, numPhBands, numWannier, numWannier}),
        KOKKOS_LAMBDA(int irP, int nu, int iw1, int iw2) {
          Kokkos::complex<double> tmp(0.0);
          for (int irE = 0; irE < numElBravaisVectors; irE++) {
            // important note: the first index iw2 runs over the k+q transform
            // while iw1 runs over k
            tmp += couplingWannier_k(irE, irP, nu, iw1, iw2) * phases_k(irE);
          }
          g1(irP, nu, iw1, iw2) = tmp;
        });
#else
#pragma omp parallel default(none) shared(g1, numPhBravaisVectors, numPhBands, numWannier, numElBravaisVectors,phases_k, couplingWannier_k)
      {
#pragma omp for collapse(4)
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          for (int nu = 0; nu < numPhBands; nu++) {
            for (int iw1 = 0; iw1 < numWannier; iw1++) {
              for (int iw2 = 0; iw2 < numWannier; iw2++) {
                g1(irP, nu, iw1, iw2) = 0.0;
              }
            }
          }
        }
#pragma omp barrier
        for (int irE = 0; irE < numElBravaisVectors; irE++) {
#pragma omp for collapse(4)
          for (int irP = 0; irP < numPhBravaisVectors; irP++) {
            for (int nu = 0; nu < numPhBands; nu++) {
              for (int iw1 = 0; iw1 < numWannier; iw1++) {
                for (int iw2 = 0; iw2 < numWannier; iw2++) {
                  g1(irP, nu, iw1, iw2) += couplingWannier_k(irE, irP, nu, iw1, iw2) * phases_k(irE);
                }
              }
            }
          }
#pragma omp barrier
        }
      }
#endif

      ComplexView4D poolElPhCached(Kokkos::ViewAllocateWithoutInitializing("poolElPhCached"), numPhBravaisVectors, numPhBands, poolNb1,
                 numWannier);

      Kokkos::parallel_for(
          "elPhCached",
          Range4D({0, 0, 0, 0},
                  {numPhBravaisVectors, numPhBands, poolNb1, numWannier}),
          KOKKOS_LAMBDA(int irP, int nu, int ib1, int iw2) {
            Kokkos::complex<double> tmp(0.0);
            for (int iw1 = 0; iw1 < numWannier; iw1++) {
              tmp += g1(irP, nu, iw1, iw2) * eigvec1_k(ib1, iw1);
            }
            poolElPhCached(irP, nu, ib1, iw2) = tmp;
          });

      if (mpi->getSize(mpi->intraPoolComm)>1) {
        // copy from accelerator to CPU
        auto poolElPhCached_h = Kokkos::create_mirror_view(poolElPhCached);
        Kokkos::deep_copy(poolElPhCached_h, poolElPhCached);

        // do an mpi->allreduce across the pool
        Eigen::Tensor<std::complex<double>, 4> reduc(numPhBravaisVectors, numPhBands, poolNb1, numWannier);
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          for (int nu = 0; nu < numPhBands; nu++) {
            for (int ib1 = 0; ib1 < poolNb1; ib1++) {
              for (int iw2 = 0; iw2 < numWannier; iw2++) {
                reduc(irP, nu, ib1, iw2) = poolElPhCached_h(irP, nu, ib1, iw2);
              }
            }
          }
        }
        mpi->allReduceSum(&reduc, mpi->intraPoolComm);
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          for (int nu = 0; nu < numPhBands; nu++) {
            for (int ib1 = 0; ib1 < poolNb1; ib1++) {
              for (int iw2 = 0; iw2 < numWannier; iw2++) {
                poolElPhCached_h(irP, nu, ib1, iw2) = reduc(irP, nu, ib1, iw2);
              }
            }
          }
        }

        // if the process owns this kpoint, copy back from CPU to accelerator
        if (mpi->getRank(mpi->intraPoolComm) == iPool) {
          Kokkos::realloc(elPhCached, numPhBravaisVectors, numPhBands, poolNb1,
                          numWannier);
          Kokkos::deep_copy(elPhCached, poolElPhCached_h);
        }
      } else { // only one element in the pool -> no need for MPI comms
        Kokkos::realloc(elPhCached, numPhBravaisVectors, numPhBands, poolNb1,
                        numWannier);
        Kokkos::deep_copy(elPhCached, poolElPhCached);

      }
    }
  }

  if (eigvec1.squaredNorm()<1.0e-8) {
    return;
  }

  // get nb2 for each ik and find the max
  // since loops and views must be rectangular, not ragged
  IntView1D nb2s_k("nb2s", numLoops);
  int nb2max = 0;
  auto nb2s_h = Kokkos::create_mirror_view(nb2s_k);
  for (int ik = 0; ik < numLoops; ik++) {
    nb2s_h(ik) = eigvecs2[ik].cols();
    if (nb2s_h(ik) > nb2max) {
      nb2max = nb2s_h(ik);
    }
  }
  Kokkos::deep_copy(nb2s_k, nb2s_h);

  IntView1D usePolarCorrections("usePolarCorrections", numLoops);
  ComplexView4D polarCorrections(Kokkos::ViewAllocateWithoutInitializing("polarCorrections"),
      numLoops, numPhBands, nb1, nb2max);
  auto usePolarCorrections_h = Kokkos::create_mirror_view(usePolarCorrections);
  auto polarCorrections_h = Kokkos::create_mirror_view(polarCorrections);

// precompute all needed polar corrections
#pragma omp parallel for default(none) shared(polarCorrections_h, nb1, usePolarCorrections_h, numLoops, q3Cs, eigvecs2, eigvecs3, usePolarCorrection, nb2s_h, numPhBands, eigvec1)
  for (int ik = 0; ik < numLoops; ik++) {
    Eigen::Vector3d q3C = q3Cs[ik];

    Eigen::MatrixXcd eigvec2 = eigvecs2[ik];
    Eigen::MatrixXcd eigvec3 = eigvecs3[ik];

    usePolarCorrections_h(ik) = usePolarCorrection && q3C.norm() > 1.0e-8;
    if (usePolarCorrections_h(ik)) {
      Eigen::Tensor<std::complex<double>, 3> polarCorrection =
          getPolarCorrection(q3C, eigvec1, eigvec2, eigvec3);
      for (int nu = 0; nu < numPhBands; nu++) {
        for (int ib1 = 0; ib1 < nb1; ib1++) {
          for (int ib2 = 0; ib2 < nb2s_h(ik); ib2++) {
            polarCorrections_h(ik, nu, ib1, ib2) =
                polarCorrection(ib2, ib1, nu);
          }
        }
      }
    }
  }
  Kokkos::deep_copy(polarCorrections, polarCorrections_h);
  Kokkos::deep_copy(usePolarCorrections, usePolarCorrections_h);

  // copy eigenvectors etc. to device
  DoubleView1D phBravaisVectorsDegeneracies_k("degens", numPhBravaisVectors);
  DoubleView2D phBravaisVectors_k("bravais", numPhBravaisVectors, 3),
      q3Cs_k("q3", numLoops, 3);
  ComplexView3D eigvecs2Dagger_k("ev2Dagger", numLoops, numWannier, nb2max),
      eigvecs3_k("ev3", numLoops, numPhBands, numPhBands);
  {
    auto eigvecs2Dagger_h = Kokkos::create_mirror_view(eigvecs2Dagger_k);
    auto eigvecs3_h = Kokkos::create_mirror_view(eigvecs3_k);
    auto q3Cs_h = Kokkos::create_mirror_view(q3Cs_k);
    auto phBravaisVectorsDegeneracies_h =
        Kokkos::create_mirror_view(phBravaisVectorsDegeneracies_k);
    auto phBravaisVectors_h = Kokkos::create_mirror_view(phBravaisVectors_k);

#pragma omp parallel for default(none) shared(eigvecs3_h, eigvecs2Dagger_h, nb2s_h, q3Cs_h, q3Cs_k, q3Cs, numLoops, numWannier, numPhBands, eigvecs2Dagger_k, eigvecs3_k, eigvecs2, eigvecs3)
    for (int ik = 0; ik < numLoops; ik++) {
      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < nb2s_h(ik); j++) {
          eigvecs2Dagger_h(ik, i, j) = std::conj(eigvecs2[ik](i, j));
        }
      }
      for (int i = 0; i < numPhBands; i++) {
        for (int j = 0; j < numPhBands; j++) {
          eigvecs3_h(ik, i, j) = eigvecs3[ik](j, i);
        }
      }
      for (int i = 0; i < numPhBands; i++) {
        for (int j = 0; j < numPhBands; j++) {
          eigvecs3_h(ik, i, j) = eigvecs3[ik](j, i);
        }
      }

      for (int i = 0; i < 3; i++) {
        q3Cs_h(ik, i) = q3Cs[ik](i);
      }
    }
    for (int i = 0; i < numPhBravaisVectors; i++) {
      phBravaisVectorsDegeneracies_h(i) = phBravaisVectorsDegeneracies(i);
      for (int j = 0; j < 3; j++) {
        phBravaisVectors_h(i, j) = phBravaisVectors(j, i);
      }
    }
    Kokkos::deep_copy(eigvecs2Dagger_k, eigvecs2Dagger_h);
    Kokkos::deep_copy(eigvecs3_k, eigvecs3_h);
    Kokkos::deep_copy(q3Cs_k, q3Cs_h);
    Kokkos::deep_copy(phBravaisVectorsDegeneracies_k,
                      phBravaisVectorsDegeneracies_h);
    Kokkos::deep_copy(phBravaisVectors_k, phBravaisVectors_h);
  }

  ComplexView2D phases("phases", numLoops, numPhBravaisVectors);
  ComplexView4D g3(Kokkos::ViewAllocateWithoutInitializing("g3"), numLoops, numPhBands, nb1, numWannier),
      g4(Kokkos::ViewAllocateWithoutInitializing("g4"), numLoops, numPhBands, nb1, numWannier),
      gFinal(Kokkos::ViewAllocateWithoutInitializing("gfinal"), numLoops, numPhBands, nb1, nb2max);
  DoubleView4D coupling_k(Kokkos::ViewAllocateWithoutInitializing("coupling"), numLoops, numPhBands, nb2max, nb1);

  Kokkos::parallel_for(
      "phases", Range2D({0, 0}, {numLoops, numPhBravaisVectors}),
      KOKKOS_LAMBDA(int ik, int irP) {
        double arg = 0.0;
        for (int j = 0; j < 3; j++) {
          arg += q3Cs_k(ik, j) * phBravaisVectors_k(irP, j);
        }
        phases(ik, irP) =
            exp(complexI * arg) / phBravaisVectorsDegeneracies_k(irP);
      });
  Kokkos::parallel_for(
      "g3", Range4D({0, 0, 0, 0}, {numLoops, numPhBands, nb1, numWannier}),
      KOKKOS_LAMBDA(int ik, int nu, int ib1, int iw2) {
        Kokkos::complex<double> tmp(0.0);
        for (int irP = 0; irP < numPhBravaisVectors; irP++) {
          tmp += phases(ik, irP) * elPhCached(irP, nu, ib1, iw2);
        }
        g3(ik, nu, ib1, iw2) = tmp;
      });

  Kokkos::parallel_for(
      "g4", Range4D({0, 0, 0, 0}, {numLoops, numPhBands, nb1, numWannier}),
      KOKKOS_LAMBDA(int ik, int nu2, int ib1, int iw2) {
        Kokkos::complex<double> tmp(0.0);
        for (int nu = 0; nu < numPhBands; nu++) {
          tmp += g3(ik, nu, ib1, iw2) * eigvecs3_k(ik, nu2, nu);
        }
        g4(ik, nu2, ib1, iw2) = tmp;
      });

  Kokkos::parallel_for(
      "gFinal", Range4D({0, 0, 0, 0}, {numLoops, numPhBands, nb1, nb2max}),
      KOKKOS_LAMBDA(int ik, int nu, int ib1, int ib2) {
        Kokkos::complex<double> tmp(0.0);
        for (int iw2 = 0; iw2 < numWannier; iw2++) {
          tmp += eigvecs2Dagger_k(ik, iw2, ib2) * g4(ik, nu, ib1, iw2);
        }
        gFinal(ik, nu, ib1, ib2) = tmp;
      });
  if (usePolarCorrection) {
    Kokkos::parallel_for(
        "correction",
        Range4D({0, 0, 0, 0}, {numLoops, numPhBands, nb1, nb2max}),
        KOKKOS_LAMBDA(int ik, int nu, int ib1, int ib2) {
          gFinal(ik, nu, ib1, ib2) += polarCorrections(ik, nu, ib1, ib2);
        });
  }
  Kokkos::parallel_for(
      "coupling", Range4D({0, 0, 0, 0}, {numLoops, numPhBands, nb2max, nb1}),
      KOKKOS_LAMBDA(int ik, int nu, int ib2, int ib1) {
        // notice the flip of 1 and 2 indices is intentional
        // coupling is |<k+q,ib2 | dV_nu | k,ib1>|^2
        auto tmp = gFinal(ik, nu, ib1, ib2);
        coupling_k(ik, nu, ib2, ib1) =
            tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
      });
  auto coupling_h = Kokkos::create_mirror_view(coupling_k);
  Kokkos::deep_copy(coupling_h, coupling_k);
#pragma omp parallel for default(none) shared(numLoops, cacheCoupling, coupling_h, nb1, nb2s_h, numPhBands)
  for (int ik = 0; ik < numLoops; ik++) {
    Eigen::Tensor<double, 3> coupling(nb1, nb2s_h(ik), numPhBands);
    for (int nu = 0; nu < numPhBands; nu++) {
      for (int ib2 = 0; ib2 < nb2s_h(ik); ib2++) {
        for (int ib1 = 0; ib1 < nb1; ib1++) {
          coupling(ib1, ib2, nu) = coupling_h(ik, nu, ib2, ib1);
        }
      }
    }
    cacheCoupling[ik] = coupling;
  }
}

Eigen::VectorXi InteractionElPhWan::getCouplingDimensions() {
  auto x = couplingWannier.dimensions();
  Eigen::VectorXi xx(5);
  for (int i : {0,1,2,3,4}) {
    xx(i) = int(x[i]);
  }
  return xx;
}
