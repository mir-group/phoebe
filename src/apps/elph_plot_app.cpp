#include "elph_plot_app.h"
#include "bandstructure.h"
#include "context.h"
#include "el_scattering.h"
#include "exceptions.h"
#include "io.h"
#include "points.h"
#include "parser.h"

#ifdef HDF5_AVAIL
#include <highfive/H5Easy.hpp>
#endif

void ElPhCouplingPlotApp::run(Context &context) {

  // load ph files
  auto t2 = Parser::parsePhHarmonic(context);
  auto crystal = std::get<0>(t2);
  auto phononH0 = std::get<1>(t2);

  // load electronic files
  auto t1 = Parser::parseElHarmonicWannier(context, &crystal);
  auto crystalEl = std::get<0>(t1);
  auto electronH0 = std::get<1>(t1);

  // load the el-ph coupling
  // Note: this file contains the number of electrons
  // which is needed to understand where to place the fermi level
  auto couplingElPh = InteractionElPhWan::parse(context, crystal, &phononH0);

  Points points(crystal);
  // decide what kind of points path we're going to use
  if (context.getG2MeshStyle() == "pointsPath") {
    points = Points(crystal, context.getPathExtrema(), context.getDeltaPath());
  }
  else { //(context.getG2MeshStyle() == "pointsMesh") { // pointsMesh is default
    points = Points(crystal, context.getKMesh());
  }
  //else {
   // Error("Select a valid g2MeshStyle "
   //     "-- options are either 'pointsMesh' or 'pointsPath'.");
  //}

  // loop over points and set up points pairs
  std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> pointsPath;
  for (int ik = 0; ik < points.getNumPoints(); ik++) {

    auto thisPoint = points.getPointCoordinates(ik, Points::cartesianCoordinates);
    std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair;

    // create a list of (k,q) pairs, where k is on a path and q is fixed
    if (context.getG2PlotStyle() == "qFixed") {
      thisPair.first = thisPoint;
      thisPair.second = context.getG2PlotFixedPoint();
    }
    // create a list of (k,q) pairs, where k is fixed and q is on the path
    else if (context.getG2PlotStyle() == "kFixed") {
      thisPair.first = context.getG2PlotFixedPoint();
      thisPair.second = thisPoint;
    }
    else { // if neither point is fixed, it's all to all
      for (int iq = 0; iq < points.getNumPoints(); iq++) {
        thisPair.first = thisPoint;
        thisPair.second = points.getPointCoordinates(iq, Points::cartesianCoordinates);;
      }
    }
    pointsPath.push_back(thisPair);
  }

  // band ranges to calculate the coupling for
  std::pair<int, int> g2PlotEl1Bands = context.getG2PlotEl1Bands();
  std::pair<int, int> g2PlotEl2Bands = context.getG2PlotEl2Bands();
  std::pair<int, int> g2PlotPhBands = context.getG2PlotPhBands();

  // Compute the coupling
  std::vector<std::complex<double>> allGs;

  // distribute over k,q pairs
  int numPairs = pointsPath.size();
  auto pairParallelIter = mpi->divideWorkIter(numPairs);
  // we calculate the coupling for each pair, flatten it, and append
  // it to allGs. Then at the end, we write this chunk to HDF5.
  for (auto iPair : pairParallelIter) {

    std::pair<Eigen::Vector3d, Eigen::Vector3d> thisPair = pointsPath[iPair];

    Eigen::Vector3d k1C = thisPair.first;
    Eigen::Vector3d q3C = thisPair.second;
    Eigen::Vector3d k2C = k1C + q3C;

    // need to get the eigenvectors at these three wavevectors
    auto t3 = electronH0.diagonalizeFromCoordinates(k1C);
    auto eigenVector1 = std::get<1>(t3);

    // second electron eigenvector
    auto t4 = electronH0.diagonalizeFromCoordinates(k2C);
    auto eigenVector2 = std::get<1>(t4);

    std::vector<Eigen::MatrixXcd> eigenVectors2;
    eigenVectors2.push_back(eigenVector2);
    std::vector<Eigen::Vector3d> k2Cs;
    k2Cs.push_back(k2C);

    // phonon eigenvectors
    auto t5 = phononH0.diagonalizeFromCoordinates(q3C);
    auto eigenVector3 = std::get<1>(t5);

    std::vector<Eigen::MatrixXcd> eigenVectors3;
    eigenVectors3.push_back(eigenVector3);
    std::vector<Eigen::Vector3d> q3Cs;
    q3Cs.push_back(q3C);

    // calculate polar correction
    std::vector<Eigen::VectorXcd> polarData;
    Eigen::VectorXcd polar = couplingElPh.polarCorrectionPart1(q3C, eigenVector3);
    polarData.push_back(polar);

    // calculate the elph coupling
    couplingElPh.cacheElPh(eigenVector1, k1C);
    couplingElPh.calcCouplingSquared(eigenVector1, eigenVectors2, eigenVectors3,
                                     q3Cs, polarData);
    auto coupling = couplingElPh.getCouplingSquared(0);

    // the coupling object is coupling at a given set of k,q, for a range of bands
    for (int ib1 = g2PlotEl1Bands.first; ib1 <= g2PlotEl1Bands.second; ib1++) {
      for (int ib2 = g2PlotEl2Bands.first; ib2 <= g2PlotEl2Bands.second; ib2++) {
        for (int ib3 = g2PlotPhBands.first; ib3 <= g2PlotPhBands.second; ib3++) {
          allGs.push_back(coupling(ib1, ib2, ib3));
        }
      }
    }
  } // close pairs loop

  // now that we've collected all the G values, we want to write them to file.

  #if defined(MPI_AVAIL) && !defined(HDF5_SERIAL)
  {
    // open the hdf5 file
    HighFive::FileAccessProps fapl;
    fapl.add(HighFive::MPIOFileAccess<MPI_Comm, MPI_Info>(MPI_COMM_WORLD, MPI_INFO_NULL));
    HighFive::File file(context.getQuantumEspressoPrefix()+".gmatrix.phoebe.hdf5",
        HighFive::File::Overwrite, fapl);

    // product of nbands1 * nbands2 * nmodes
    size_t bandProd = (g2PlotEl1Bands.second - g2PlotEl1Bands.first) *
                (g2PlotEl2Bands.second - g2PlotEl2Bands.first) *
                (g2PlotPhBands.second - g2PlotPhBands.first);

    unsigned int globalSize = numPairs * bandProd;

    // Create the data-space to write g to
    std::vector<size_t> dims(2);
    dims[0] = 1;
    dims[1] = size_t(globalSize);
    HighFive::DataSet dgmat = file.createDataSet<std::complex<double>>(
          "/gMat", HighFive::DataSpace(dims));

    // start point and the number of the total number of elements
    // to be written by this process
    size_t start = mpi->divideWorkIter(numPairs)[0] * bandProd;
    size_t offset = start;

    // Note: HDF5 < v1.10.2 cannot write datasets larger than 2 Gbs
    // ( due to max(int 32 bit))/1024^3 = 2Gb overflowing in MPI)
    // In order to be compatible with older versions, we split the tensor
    // into smaller chunks and write them to separate datasets
    // slower, but it will work more often.

    // maxSize represents 2GB worth of std::complex<doubles>, since that's
    // what we write
    auto maxSize = int(pow(1000, 3)) / sizeof(std::complex<double>);
    size_t smallestSize = bandProd; // 1 point pair
    std::vector<int> bunchSizes;

    // determine the size of each bunch of electronic bravais vectors
    // the BunchSizes vector tells us how many are in each set
    int numPairsBunch = mpi->divideWorkIter(numPairs).back() + 1 - mpi->divideWorkIter(numPairs)[0];

    int bunchSize = 0;
    for (int i = 0; i < numPairsBunch; i++) {
      bunchSize++;
      // this bunch is as big as possible, stop adding to it
      if ((bunchSize + 1) * smallestSize > maxSize) {
        bunchSizes.push_back(bunchSize);
        bunchSize = 0;
      }
    }
    // push the last one, no matter the size, to the list of bunch sizes
    bunchSizes.push_back(bunchSize);

    // determine the number of bunches. The last bunch may be smaller
    // than the rest
    int numDatasets = bunchSizes.size();

    // we now loop over these data sets, and write each chunk of
    // bravais vectors in parallel
    int netOffset = 0;// offset from first bunch in this set to current bunch
    for (int iBunch = 0; iBunch < numDatasets; iBunch++) {

      // we need to determine the start, stop and offset of this
      // sub-slice of the dataset available to this process
      size_t bunchElements = bunchSizes[iBunch] * smallestSize;
      //size_t bunchStart = start + netOffset;
      size_t bunchOffset = offset + netOffset;
      netOffset += bunchElements;

      // Each process writes to hdf5
      // The format is ((startRow,startCol),(numRows,numCols)).write(data)
      // Because it's a vector (1 row) all processes write to row=0, col=startPoint
      // with nRows = 1, nCols = number of items this process will write.
      dgmat.select({0, bunchOffset}, {1, bunchElements}).write_raw(&allGs);
    }
  }
  #else
  Error("You cannot output the elph matrix elements to HDF5 because your copy of \n"
        "Phoebe has not been compiled with HDF5 support,\n"
        "or has been compiled with serial HDF5.");
  #endif
}

void ElPhCouplingPlotApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  throwErrorIfUnset(context.getPhFC2FileName(), "phFC2FileName");
  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getElphFileName(), "elphFileName");
}
