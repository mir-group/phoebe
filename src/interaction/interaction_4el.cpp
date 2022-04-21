#include "interaction_4el.h"
#include <Kokkos_Core.hpp>
#include <fstream>

#ifdef HDF5_AVAIL
#include <Kokkos_ScatterView.hpp>
#include <highfive/H5Easy.hpp>
#endif

// default constructor
Interaction4El::Interaction4El(
    Crystal &crystal_,
    const Eigen::Tensor<std::complex<double>, 7> &couplingWannier_,
    const Eigen::MatrixXd &elBravaisVectors,
    const Eigen::VectorXd &elBravaisVectorsDegeneracies)
    : crystal(crystal_) {

  int numElBravaisVectors = couplingWannier_.dimension(0);
  int numWannier = couplingWannier_.dimension(4);

  HostComplexView7D couplingWannier_h((Kokkos::complex<double>*) couplingWannier_.data(),
                                      numElBravaisVectors, numElBravaisVectors, numElBravaisVectors,
                                      numWannier, numWannier, numWannier, numWannier);
  couplingWannier_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), couplingWannier_h);
  Kokkos::deep_copy(couplingWannier_d, couplingWannier_h);

  HostDoubleView2D elBravaisVectors_h((double*) elBravaisVectors.data(),
                                      numElBravaisVectors, 3);
  elBravaisVectors_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), elBravaisVectors_h);
  Kokkos::deep_copy(elBravaisVectors_d, elBravaisVectors_h);

  HostDoubleView1D elBravaisVectorsDegeneracies_h((double*) elBravaisVectorsDegeneracies.data(),
                                      numElBravaisVectors);
  elBravaisVectorsDegeneracies_d = Kokkos::create_mirror_view(Kokkos::DefaultExecutionSpace(), elBravaisVectorsDegeneracies_h);
  Kokkos::deep_copy(elBravaisVectorsDegeneracies_d, elBravaisVectorsDegeneracies_h);
}

// copy constructor
Interaction4El::Interaction4El(const Interaction4El &that)
    : crystal(that.crystal),
      numWannier(that.numWannier),
      numElBravaisVectors(that.numElBravaisVectors),
      cacheCoupling(that.cacheCoupling),
      elPhCached1(that.elPhCached1),
      elPhCached2(that.elPhCached2),
      couplingWannier_d(that.couplingWannier_d),
      elBravaisVectors_d(that.elBravaisVectors_d),
      elBravaisVectorsDegeneracies_d(that.elBravaisVectorsDegeneracies_d) {}

// assignment operator
Interaction4El &Interaction4El::operator=(const Interaction4El &that) {
  if (this != &that) {
    crystal = that.crystal;
    numWannier = that.numWannier;
    numElBravaisVectors = that.numElBravaisVectors;
    cacheCoupling = that.cacheCoupling;
    elPhCached1 = that.elPhCached1;
    elPhCached2 = that.elPhCached2;
    couplingWannier_d = that.couplingWannier_d;
    elBravaisVectors_d = that.elBravaisVectors_d;
    elBravaisVectorsDegeneracies_d = that.elBravaisVectorsDegeneracies_d;
  }
  return *this;
}

Eigen::Tensor<double, 4> Interaction4El::getCouplingSquared(const int &ik3) {
  return cacheCoupling[ik3];
}

void Interaction4El::calcCouplingSquared(
    const std::vector<Eigen::MatrixXcd> &eigvecs3,
    const std::vector<Eigen::MatrixXcd> &eigvecs4,
    const std::vector<Eigen::Vector3d> &k3Cs,
    const std::vector<Eigen::Vector3d> &k4Cs) {

  auto numLoops = int(eigvecs3.size());

  auto elPhCached2 = this->elPhCached2;
  int numElBravaisVectors = this->numElBravaisVectors;
  auto elBravaisVectors_d = this->elBravaisVectors_d;
  auto elBravaisVectorsDegeneracies_d = this->elBravaisVectorsDegeneracies_d;

  int nb1 = elPhCached2.extent(1);
  int nb2 = elPhCached2.extent(2);

  // get nb2 for each ik and find the max
  // since loops and views must be rectangular, not ragged
  IntView1D nb3s_d("nb3s", numLoops);
  IntView1D nb4s_d("nb4s", numLoops);
  int nb3max = 0;
  int nb4max = 0;
  auto nb3s_h = Kokkos::create_mirror_view(nb3s_d);
  auto nb4s_h = Kokkos::create_mirror_view(nb4s_d);
  for (int ik = 0; ik < numLoops; ik++) {
    nb3s_h(ik) = int(eigvecs3[ik].cols());
    nb4s_h(ik) = int(eigvecs4[ik].cols());
    if (nb3s_h(ik) > nb3max) {
      nb3max = nb3s_h(ik);
    }
    if (nb4s_h(ik) > nb4max) {
      nb4max = nb4s_h(ik);
    }
  }
  Kokkos::deep_copy(nb3s_d, nb3s_h);
  Kokkos::deep_copy(nb4s_d, nb4s_h);

  // copy eigenvectors etc. to device
  DoubleView2D k3Cs_d("k3", numLoops, 3);
  ComplexView3D eigvecs3_d("ev3", numLoops, numWannier, nb3max);
  ComplexView3D eigvecs4_d("ev4", numLoops, numWannier, nb4max);
  {
    auto eigvecs3_h = Kokkos::create_mirror_view(eigvecs3_d);
    auto eigvecs4_h = Kokkos::create_mirror_view(eigvecs4_d);
    auto k3Cs_h = Kokkos::create_mirror_view(k3Cs_d);

#pragma omp parallel for
    for (int ik = 0; ik < numLoops; ik++) {

      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < nb3s_h(ik); j++) {
          eigvecs3_h(ik, i, j) = eigvecs3[ik](i, j);
        }
      }
      for (int i = 0; i < numWannier; i++) {
        for (int j = 0; j < nb4s_h(ik); j++) {
          eigvecs4_h(ik, i, j) = eigvecs4[ik](i, j);
        }
      }

      for (int i = 0; i < 3; i++) {
        k3Cs_h(ik, i) = k3Cs[ik](i);
      }
    }
    Kokkos::deep_copy(eigvecs3_d, eigvecs3_h);
    Kokkos::deep_copy(eigvecs4_d, eigvecs4_h);
    Kokkos::deep_copy(k3Cs_d, k3Cs_h);
  }

  // now we finish the Wannier transform. We have to do the Fourier transform
  // on the lattice degrees of freedom, and then do two rotations (at k2 and q)
  ComplexView2D phases("phases", numLoops, numElBravaisVectors);
  Kokkos::complex<double> complexI(0.0, 1.0);
  Kokkos::parallel_for(
      "phases", Range2D({0, 0}, {numLoops, numElBravaisVectors}),
      KOKKOS_LAMBDA(int ik, int irE3) {
        double arg = 0.0;
        for (int j = 0; j < 3; j++) {
          arg += k3Cs_d(ik, j) * elBravaisVectors_d(irE3, j);
        }
        phases(ik, irE3) =
            exp(complexI * arg) / elBravaisVectorsDegeneracies_d(irE3);
      });

  // Last Fourier transform term
  ComplexView5D g3(Kokkos::ViewAllocateWithoutInitializing("g3"), numLoops, nb1, nb2, numWannier, numWannier);
  Kokkos::parallel_for(
      "g3", Range5D({0, 0, 0, 0, 0}, {numLoops, nb1, nb2, numWannier, numWannier}),
      KOKKOS_LAMBDA(int ik, int ib1, int ib2, int iw3, int iw4) {
        Kokkos::complex<double> tmp(0., 0.);
        for (int irE3 = 0; irE3 < numElBravaisVectors; irE3++) {
          tmp += phases(ik, irE3) * elPhCached2(irE3, ib1, ib2, iw3, iw4);
        }
        g3(ik, ib1, ib2, iw3, iw4) = tmp;
      });
  Kokkos::realloc(phases, 0, 0);

  // now we rotate with the eigenvectors
  ComplexView5D g4(Kokkos::ViewAllocateWithoutInitializing("g4"), numLoops, nb1, nb2, nb3max, numWannier);
  Kokkos::parallel_for(
      "g4", Range5D({0, 0, 0, 0, 0}, {numLoops, nb1, nb2, nb3max, numWannier}),
      KOKKOS_LAMBDA(int ik, int ib1, int ib2, int ib3, int iw4) {
        Kokkos::complex<double> tmp(0.,0.);
        for (int iw3 = 0; iw3 < numWannier; iw3++) {
          tmp += g3(ik, ib1, ib2, iw3, iw4) * eigvecs3_d(ik, ib3, iw3);
        }
        g4(ik, ib1, ib2, ib3, iw4) = tmp;
      });
  Kokkos::realloc(g3, 0, 0, 0, 0, 0);

  ComplexView5D gFinal(Kokkos::ViewAllocateWithoutInitializing("gFinal"), numLoops, nb1, nb2, nb3max, nb4max);
  Kokkos::parallel_for(
      "gFinal", Range5D({0, 0, 0, 0, 0}, {numLoops, nb1, nb2, nb3max, nb4max}),
      KOKKOS_LAMBDA(int ik, int ib1, int ib2, int ib3, int ib4) {
        Kokkos::complex<double> tmp(0.,0.);
        for (int iw4 = 0; iw4 < numWannier; iw4++) {
          tmp += g4(ik, ib1, ib2, ib3, iw4) * eigvecs4_d(ik, ib4, iw4);
        }
        gFinal(ik, ib1, ib2, ib3, ib4) = tmp;
      });
  Kokkos::realloc(g4, 0, 0, 0, 0, 0);

  // finally, compute |g|^2 from g
  DoubleView5D coupling_d(Kokkos::ViewAllocateWithoutInitializing("coupling"), numLoops, nb1, nb2, nb3max, nb4max);
  Kokkos::parallel_for(
      "coupling", Range5D({0, 0, 0, 0, 0}, {numLoops, nb1, nb2, nb3max, nb4max}),
      KOKKOS_LAMBDA(int ik, int ib1, int ib2, int ib3, int ib4) {
        // notice the flip of 1 and 2 indices is intentional
        // coupling is |<k+q,ib2 | dV_nu | k,ib1>|^2
        auto tmp = gFinal(ik, ib1, ib2, ib3, ib4);
        coupling_d(ik, ib1, ib2, ib3, ib4) =
            tmp.real() * tmp.real() + tmp.imag() * tmp.imag();
      });
  Kokkos::realloc(gFinal, 0, 0, 0, 0, 0);

  // now, copy results back to the CPU
  cacheCoupling.resize(0);
  cacheCoupling.resize(numLoops);
  auto coupling_h = Kokkos::create_mirror_view(coupling_d);
  Kokkos::deep_copy(coupling_h, coupling_d);
#pragma omp parallel for
  for (int ik = 0; ik < numLoops; ik++) {
    Eigen::Tensor<double, 4> coupling(nb1, nb2, nb3s_h(ik), nb4s_h(ik));
    for (int ib4 = 0; ib4 < nb4s_h(ik); ib4++) {
      for (int ib3 = 0; ib3 < nb3s_h(ik); ib3++) {
        for (int ib2 = 0; ib2 < nb2; ib2++) {
          for (int ib1 = 0; ib1 < nb1; ib1++) {
            coupling(ib1, ib2, ib3, ib4) = coupling_h(ik, ib1, ib2, ib3, ib4);
          }
        }
      }
    }
    // and we save the coupling |g|^2 it for later
    cacheCoupling[ik] = coupling;
  }
}

void Interaction4El::cache1stEl(const Eigen::MatrixXcd &eigvec1, const Eigen::Vector3d &k1C) {
  Kokkos::complex<double> complexI(0.0, 1.0);

  // note: when Kokkos is compiled with GPU support, we must create elPhCached
  // and other variables as local, so that Kokkos correctly allocates these
  // quantities on the GPU. At the end of this function, elPhCached must be
  // 'copied' back into this->elPhCached. Note that no copy actually is done,
  // since Kokkos::View works similarly to a shared_ptr.
  auto elPhCached1 = this->elPhCached1;
  int numElBravaisVectors = this->numElBravaisVectors;
  int numWannier = this->numWannier;
  auto elBravaisVectors_d = this->elBravaisVectors_d;
  auto elBravaisVectorsDegeneracies_d = this->elBravaisVectorsDegeneracies_d;

  // move eigenvectors and wavevectors to device
  int nb1 = int(eigvec1.cols());
  ComplexView2D eigvec1_d("ev1", nb1, numWannier);
  DoubleView1D k1C_d("k", 3);
  {
    HostComplexView2D eigvec1_h((Kokkos::complex<double>*) eigvec1.data(), nb1, numWannier);
    HostDoubleView1D k1C_h((double*) k1C.data(), 3);
    Kokkos::deep_copy(eigvec1_d, eigvec1_h);
    Kokkos::deep_copy(k1C_d, k1C_h);
  }

  // first we precompute the phases
  ComplexView1D phases_k("phases", numElBravaisVectors);
  Kokkos::parallel_for(
      "phases_k", numElBravaisVectors,
      KOKKOS_LAMBDA(int irE) {
        double arg = 0.0;
        for (int j = 0; j < 3; j++) {
          arg += k1C_d(j) * elBravaisVectors_d(irE, j);
        }
        phases_k(irE) =
            exp(complexI * arg) / elBravaisVectorsDegeneracies_d(irE);
      });

  // fourier transform on 1st coordinate
  ComplexView6D preCache1("preCache1", numElBravaisVectors, numElBravaisVectors,
                          numWannier, numWannier, numWannier, numWannier);
  Kokkos::parallel_for(
      "preCache1",
      Range6D({0, 0, 0, 0, 0, 0},
              {numElBravaisVectors, numElBravaisVectors, numWannier, numWannier, numWannier, numWannier}),
      KOKKOS_LAMBDA(int irE2, int irE3, int iw1, int iw2, int iw3, int iw4) {
        Kokkos::complex<double> tmp(0.0);
        for (int irE1 = 0; irE1 < numElBravaisVectors; irE1++) {
          // important note: the first index iw2 runs over the k+q transform
          // while iw1 runs over k
          tmp += couplingWannier_d(irE1, irE2, irE3, iw1, iw2, iw3, iw4) * phases_k(irE1);
        }
        preCache1(irE2, irE3, iw1, iw2, iw3, iw4) = tmp;
      });

  // rotate eigenvector
  Kokkos::realloc(elPhCached1, numElBravaisVectors, numElBravaisVectors, nb1,
                  numWannier, numWannier, numWannier);
  Kokkos::parallel_for(
      "cache1",
      Range6D({0, 0, 0, 0, 0, 0},
              {numElBravaisVectors, numElBravaisVectors, nb1, numWannier, numWannier, numWannier}),
      KOKKOS_LAMBDA(int irE2, int irE3, int ib1, int iw2, int iw3, int iw4) {
        Kokkos::complex<double> tmp(0.0);
        for (int iw1 = 0; iw1 < numWannier; iw1++) {
          tmp += preCache1(irE2, irE3, iw1, iw2, iw3, iw4) * eigvec1_d(ib1, iw1);
        }
        elPhCached1(irE2, irE3, ib1, iw2, iw3, iw4) = tmp;
      });
  this->elPhCached1 = elPhCached1;
}

void Interaction4El::cache2ndEl(const Eigen::MatrixXcd &eigvec2, const Eigen::Vector3d &k2C) {
  Kokkos::complex<double> complexI(0.0, 1.0);

  // note: when Kokkos is compiled with GPU support, we must create elPhCached
  // and other variables as local, so that Kokkos correctly allocates these
  // quantities on the GPU. At the end of this function, elPhCached must be
  // 'copied' back into this->elPhCached. Note that no copy actually is done,
  // since Kokkos::View works similarly to a shared_ptr.
  auto elPhCached1 = this->elPhCached1;
  auto elPhCached2 = this->elPhCached2;
  int numElBravaisVectors = this->numElBravaisVectors;
  int numWannier = this->numWannier;
  auto elBravaisVectors_d = this->elBravaisVectors_d;
  auto elBravaisVectorsDegeneracies_d = this->elBravaisVectorsDegeneracies_d;

  int nb1 = elPhCached1.extent(0);

  // move eigenvectors and wavevectors to device
  int nb2 = int(eigvec2.cols());
  ComplexView2D eigvec2_d("ev2", nb2, numWannier);
  DoubleView1D k2C_d("k", 3);
  {
    HostComplexView2D eigvec2_h((Kokkos::complex<double>*) eigvec2.data(), nb2, numWannier);
    HostDoubleView1D k2C_h((double*) k2C.data(), 3);
    Kokkos::deep_copy(eigvec2_d, eigvec2_h);
    Kokkos::deep_copy(k2C_d, k2C_h);
  }

  // first we precompute the phases
  ComplexView1D phases_k("phases", numElBravaisVectors);
  Kokkos::parallel_for(
      "phases_k", numElBravaisVectors,
      KOKKOS_LAMBDA(int irE) {
        double arg = 0.0;
        for (int j = 0; j < 3; j++) {
          arg += k2C_d(j) * elBravaisVectors_d(irE, j);
        }
        phases_k(irE) =
            exp(complexI * arg) / elBravaisVectorsDegeneracies_d(irE);
      });

  // fourier transform on 1st coordinate
  ComplexView5D preCache2("preCache2", numElBravaisVectors,
                          nb1, numWannier, numWannier, numWannier);
  Kokkos::parallel_for(
      "preCache2",
      Range5D({0, 0, 0, 0, 0},
              {numElBravaisVectors, numWannier, numWannier, numWannier, numWannier}),
      KOKKOS_LAMBDA(int irE3, int ib1, int iw2, int iw3, int iw4) {
        Kokkos::complex<double> tmp(0.0);
        for (int irE2 = 0; irE2 < numElBravaisVectors; irE2++) {
          // important note: the first index iw2 runs over the k+q transform
          // while iw1 runs over k
          tmp += elPhCached1(irE2, irE3, ib1, iw2, iw3, iw4) * phases_k(irE2);
        }
        preCache2(irE3, ib1, iw2, iw3, iw4) = tmp;
      });

  // rotate eigenvector
  Kokkos::realloc(elPhCached2, numElBravaisVectors, nb1, nb2,
                  numWannier, numWannier);
  Kokkos::parallel_for(
      "cache2",
      Range5D({0, 0, 0, 0, 0},
              {numElBravaisVectors, nb1, nb2, numWannier, numWannier}),
      KOKKOS_LAMBDA(int irE3, int ib1, int ib2, int iw3, int iw4) {
        Kokkos::complex<double> tmp(0.0);
        for (int iw2 = 0; iw2 < numWannier; iw2++) {
          tmp += preCache2(irE3, ib1, iw2, iw3, iw4) * eigvec2_d(ib2, iw2);
        }
        elPhCached2(irE3, ib1, ib2, iw3, iw4) = tmp;
      });
  this->elPhCached2 = elPhCached2;
}

#ifdef HDF5_AVAIL

std::tuple<int, int, Eigen::MatrixXd, Eigen::VectorXd>
    parseHeaderHDF5(Context &context) {
  std::string fileName = context.getElphFileName();

  int numWannier;
  int numElBravaisVectors;
  // suppress initialization warning
  numWannier = 0;
  numElBravaisVectors = 0;
  Eigen::MatrixXd elBravaisVectors_;
  Eigen::VectorXd elBravaisVectorsDegeneracies_;
  Eigen::Tensor<std::complex<double>, 7> couplingWannier_;
  std::vector<size_t> localElVectors;

  try {
    // Use MPI head only to read in the small data structures
    // then distribute them below this
    if (mpi->mpiHeadPool()) {
      // need to open the files differently if MPI is available or not
      // NOTE: do not remove the braces inside this if -- the file must
      // go out of scope, so that it can be reopened for parallel
      // read in the next block.
      {
        // Open the HDF5 ElPh file
        HighFive::File file(fileName, HighFive::File::ReadOnly);

        // read in the number of phonon and electron bands
        HighFive::DataSet dnWannier = file.getDataSet("/numWannier");
        dnWannier.read(numWannier);

        // read electron Bravais lattice vectors and degeneracies
        HighFive::DataSet delDegeneracies = file.getDataSet("/elDegeneracies");
        delDegeneracies.read(elBravaisVectorsDegeneracies_);
        numElBravaisVectors = int(elBravaisVectorsDegeneracies_.size());

        HighFive::DataSet delbravais = file.getDataSet("/elBravaisVectors");
        delbravais.read(elBravaisVectors_);
      }
    }
    // broadcast to all MPI processes
    mpi->bcast(&numWannier);
    mpi->bcast(&numElBravaisVectors);

    if (!mpi->mpiHeadPool()) {// head already allocated these
      elBravaisVectors_.resize(3, numElBravaisVectors);
      elBravaisVectorsDegeneracies_.resize(numElBravaisVectors);
      couplingWannier_.resize(numWannier, numWannier, numWannier, numWannier,
                              numElBravaisVectors, numElBravaisVectors,
                              numElBravaisVectors);
    }
    mpi->bcast(&elBravaisVectors_, mpi->interPoolComm);
    mpi->bcast(&elBravaisVectorsDegeneracies_, mpi->interPoolComm);
  } catch (std::exception &error) {
    Error("Issue reading elph Wannier representation from hdf5.");
  }

  return std::make_tuple(numWannier, numElBravaisVectors, elBravaisVectors_,
                         elBravaisVectorsDegeneracies_);
}

// specific parse function for the case where parallel HDF5 is available
Interaction4El parseHDF5(Context &context, Crystal &crystal) {
  std::string fileName = context.get4ElFileName();

  auto t = parseHeaderHDF5(context);
  int numWannier = std::get<0>(t);
  int totalNumElBravaisVectors = std::get<1>(t);
  Eigen::MatrixXd elBravaisVectors_ = std::get<2>(t);
  Eigen::VectorXd elBravaisVectorsDegeneracies_ = std::get<3>(t);
  int numElBravaisVectors = elBravaisVectorsDegeneracies_.size();

  Eigen::Tensor<std::complex<double>, 7> couplingWannier_;

  try {
    // Define the eph matrix element containers

    // This is broken into parts, otherwise it can overflow if done all at once
    size_t totElems = numWannier * numWannier * numWannier * numWannier;
    totElems *= numElBravaisVectors;
    totElems *= numElBravaisVectors;
    totElems *= numElBravaisVectors;

    // user info about memory
    {
      std::complex<double> cx;
      auto x = double(totElems / pow(1024., 3) * sizeof(cx));
      if (mpi->mpiHead()) {
        std::cout << "Allocating " << x
                  << " (GB) (per MPI process) for the el-ph coupling matrix."
                  << std::endl;
      }
    }

    couplingWannier_.resize(numWannier, numWannier, numWannier, numWannier,
                            numElBravaisVectors, numElBravaisVectors, numElBravaisVectors);
    couplingWannier_.setZero();

    // Set up buffer to receive full matrix data
    size_t sliceElements = pow(numWannier,3) * pow(numElBravaisVectors,2);
    Eigen::VectorXcd slice(sliceElements);

    HighFive::File file(fileName, HighFive::File::ReadOnly);
    for (int irE1 : mpi->divideWorkIter(numElBravaisVectors)) {
      std::string datasetName = "/gWannier_" + std::to_string(irE1);
      // Set up dataset for gWannier
      HighFive::DataSet dslice = file.getDataSet(datasetName);
      dslice.read(slice);

      // Map the flattened matrix back to tensor structure
      Eigen::TensorMap<Eigen::Tensor<std::complex<double>, 6>> sliceTemp(
          slice.data(), numWannier, numWannier, numWannier, numWannier, numElBravaisVectors, numElBravaisVectors);

      for (int irE2 = 0; irE2 < numElBravaisVectors; irE2++) {
        for (int irE3 = 0; irE3 < numElBravaisVectors; irE3++) {
          for (int iw1 = 0; iw1 < numWannier; iw1++) {
            for (int iw2 = 0; iw2 < numWannier; iw2++) {
              for (int iw3 = 0; iw3 < numWannier; iw3++) {
                for (int iw4 = 0; iw4 < numWannier; iw4++) {
                  couplingWannier_(iw4, iw3, iw2, iw1, irE3, irE2, irE1)
                      = sliceTemp(iw4, iw3, iw2, iw1, irE3, irE2);
                }
              }
            }
          }
        }
      }
    }
    mpi->allReduceSum(&couplingWannier_);

  } catch (std::exception &error) {
    Error("Issue reading el-el Wannier representation from hdf5.");
  }

  Interaction4El output(crystal, couplingWannier_, elBravaisVectors_,
                        elBravaisVectorsDegeneracies_);
  return output;
}

#endif

// General parse function
Interaction4El Interaction4El::parse(Context &context, Crystal &crystal) {
  if (mpi->mpiHead()) {
    std::cout << "\nStarted parsing of el-el interaction." << std::endl;
  }
#ifdef HDF5_AVAIL
  auto output = parseHDF5(context, crystal);
#else
  Error("Didn't implement 4-el coupling parsing without HDF5");
  // auto output = parseNoHDF5(context, crystal);
#endif
  if (mpi->mpiHead()) {
    std::cout << "Finished parsing of el-ph interaction." << std::endl;
  }
  return output;
}
