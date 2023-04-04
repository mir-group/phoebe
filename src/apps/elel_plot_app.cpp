#include "elel_plot_app.h"
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

void ElElCouplingPlotApp::run(Context &context) {

  if(mpi->hasPools()) {
    Error("Cannot currently run el-el coupling plot app with\n"
        "MPI pool size (-ps) greater than 1. Please run without pools\n"
        "or let the developers know you need this feature.");
  }

  // load electronic files
  auto t0 = Parser::parseElHarmonicWannier(context);
  auto crystal = std::get<0>(t0);
  auto electronH0 = std::get<1>(t0);

  // load the coupling
  Interaction4El coupling4El = Interaction4El::parse(context, crystal);

  if(mpi->mpiHead()) std::cout << "Done setting up el-el coupling." << std::endl;

  // deciding how to specify the mesh here will be weirder --
  // it's a 5pt diagram, so 3? independent indices..
  // if we fix k1 and k2...
  // we have k1, k2, exciton Q, k3, k4
  // k1 -- independent
  // k2 -- independent
  // Q -- independent
  //k3 = k1 - Q;
  //k4 = k2 + Q  OR Eigen::Vector3d k4CTemp = k1C + k2C - k3C;
  // plus b1, b2, b3, b4

/*k3-k1＝q
k1        0.125 -1.86265e-08            0
k2         0.25 -3.72529e-08            0
k3        0.375 -1.49012e-08            0
Q      -0.125 1.86265e-08           0
q = k3 - k1 = 0.25 0 0
(0.000116351,0.000201525)
 */

  // for now we set up the coupling points by hand. Later will bring in a way to read this from input
  // for now, let's do ... coupling between fixed K1, K2, and a specified uniform kmesh and do coupling between
  // set up k-points
  //
  Eigen::Vector3i kMesh = context.getKMesh();
  Points kPoints(crystal, kMesh);
  int numKPoints = kMesh.prod();
  // TODO any way to enforce this is crystal?
  Eigen::Vector3d k1 = context.getG2PlotFixedPoint();
  k1 = kPoints.crystalToCartesian(k1);
  // TODO give some way to specify this by input

  Eigen::Vector3d k2 = {0.25, 0.25, 0.0};
  k2 = kPoints.crystalToCartesian(k2);

  //Eigen::Vector3d k1 = {0.125,0.0,0.0};
  //Eigen::Vector3d k2 = {0.25, 0.0,0.0};
  //Eigen::Vector3d k3 = {0.375,0.0,0.0};
  //k1 = kPoints.crystalToCartesian(k1);
  //k2 = kPoints.crystalToCartesian(k2);
  //k3 = points.crystalToCartesian(k3);

  // loop over points and set up points triplets
  std::vector<std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d>> pointsTriplets;
  for (int ik3 = 0; ik3 < numKPoints; ik3++) {

    Eigen::Vector3d k3 = kPoints.getPointCoordinates(ik3, Points::cartesianCoordinates);
    //std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> thisTriplet;
    pointsTriplets.push_back(std::make_tuple(k1, k2, k3));

  }

  // band ranges to calculate the coupling for
  // TODO add PlotElBand ranges variables for 3 and 4? Or just use one and enforce that?
  std::pair<int, int> g2PlotEl1Bands = std::make_pair(0,1); //context.getG2PlotEl1Bands();
  std::pair<int, int> g2PlotEl2Bands = std::make_pair(0,1); //context.getG2PlotEl2Bands();
  std::pair<int, int> g2PlotEl3Bands = std::make_pair(0,1); //context.getG2PlotPhBands();
  std::pair<int, int> g2PlotEl4Bands = std::make_pair(0,1); //context.getG2PlotPhBands();

  // Compute the coupling --------------------------------------------------
  // we calculate the coupling for each pair, flatten it, and append
  // it to allMatrixElementsSq. Then at the end, we write this chunk to HDF5.
  std::vector<double> allMatrixElementsSq;

  // distribute over triplets
  int numTriplets = pointsTriplets.size();
  auto pointsParallelIter = mpi->divideWorkIter(numTriplets);

  if(mpi->mpiHead())
    std::cout << "\nCoupling requested for " << numTriplets << " k1,k2,k3 set." << std::endl;

  LoopPrint loopPrint("calculating coupling", "triplets on this process", pointsParallelIter.size());
  for (auto iTriplet : pointsParallelIter) {

    loopPrint.update();

    std::tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> thisTriplet = pointsTriplets[iTriplet];

    // set up the points as cartesian coordinate vectors
    // we set these up as yambo does (rather than how we do this for scattering)
    // Where Yambo sets up the kpoint relations:
    // k1 - k2 = k3 - k4 = Q (where Q is the exciton momentum)
    // but we are using in our notes:
    // k1 + k2 = k3 + k4
    Eigen::Vector3d k1 = std::get<0>(thisTriplet);
    Eigen::Vector3d k2 = std::get<1>(thisTriplet);
    Eigen::Vector3d k3 = std::get<2>(thisTriplet);
    Eigen::Vector3d k4 = k1 + k2 - k3;

    std::cout << "k3 " << k3.transpose() << std::endl;

    // need to get the eigenvectors at these three wavevectors
    auto t1 = electronH0.diagonalizeFromCoordinates(k1);
    auto eigenVector1 = std::get<1>(t1);

    // second electron eigenvector
    auto t2 = electronH0.diagonalizeFromCoordinates(k2);
    auto eigenVector2 = std::get<1>(t2);

    // third electron eigenvector
    auto t3 = electronH0.diagonalizeFromCoordinates(k3);
    auto eigenVector3 = std::get<1>(t3);

    // fourth electron eigenvector
    auto t4 = electronH0.diagonalizeFromCoordinates(k4);
    auto eigenVector4 = std::get<1>(t4);

    std::vector<Eigen::MatrixXcd> eigenVectors3;
    eigenVectors3.push_back(eigenVector3);
    std::vector<Eigen::Vector3d> k3s;
    k3s.push_back(k3);

    std::vector<Eigen::MatrixXcd> eigenVectors4;
    eigenVectors4.push_back(eigenVector4);
    std::vector<Eigen::Vector3d> k4s;
    k4s.push_back(k4);

    // calculate the coupling
    coupling4El.cache1stEl(eigenVector1, k1);
    coupling4El.cache2ndEl(eigenVector2, k2);
    coupling4El.calcCouplingSquared(eigenVectors3, eigenVectors4, k3s, k4s);
    // normally this function requires a kpoint index, but we only had coupling for 1
    // kpoint set, so here it's 0.
    auto coupling = coupling4El.getCouplingSquared(0);
    //
    Eigen::Tensor<double, 4>& couplingA = std::get<0>(coupling);

    // the coupling object is coupling at a given set of k,q, for a range of bands
    // band ranges are inclusive of start and finish ones
    for (int ib1 = g2PlotEl1Bands.first; ib1 <= g2PlotEl1Bands.second; ib1++) {
      for (int ib2 = g2PlotEl2Bands.first; ib2 <= g2PlotEl2Bands.second; ib2++) {
        for (int ib3 = g2PlotEl3Bands.first; ib3 <= g2PlotEl3Bands.second; ib3++) {
          for (int ib4 = g2PlotEl4Bands.first; ib4 <= g2PlotEl4Bands.second; ib4++) {
            allMatrixElementsSq.push_back(couplingA(ib1, ib2, ib3, ib4));
            //if(ib1 == 1  && ib2 == 2 && ib3 == 1) {
              //if(mpi->mpiHead()) std::cout << "coupling " << std::setprecision(8) << coupling(ib1, ib2, ib3) << " " << k1C.transpose() << std::endl;
           // }
          }
        }
      }
    }
  } // close pairs loop
  mpi->barrier();
  loopPrint.close();

  // now that we've collected all the coupling values, we want to write them to file.
  if(mpi->mpiHead())
    std::cout << "\nFinished calculating coupling, writing to file." << std::endl;

  std::string outFileName = "coupling.4el.phoebe.hdf5";
  std::remove(&outFileName[0]);

  // product of nbands1 * nbands2 * nmodes -- + 1 is because range is inclusive
  size_t bandProd = (g2PlotEl1Bands.second - g2PlotEl1Bands.first + 1) *
              (g2PlotEl2Bands.second - g2PlotEl2Bands.first + 1) *
              (g2PlotEl3Bands.second - g2PlotEl3Bands.first + 1) *
              (g2PlotEl4Bands.second - g2PlotEl4Bands.first + 1);

  #if defined(HDF5_AVAIL)
  try {
  #if defined(MPI_AVAIL) && !defined(HDF5_SERIAL)
  { // need open/close braces so that the HDF5 file goes out of scope

    // open the hdf5 file
    HighFive::File file(
          outFileName, HighFive::File::Overwrite,
          HighFive::MPIOFileDriver(MPI_COMM_WORLD, MPI_INFO_NULL));

    unsigned int globalSize = numTriplets * bandProd;

    // Create the data-space to write g to
    std::vector<size_t> dims(2);
    dims[0] = 1;
    dims[1] = size_t(globalSize);
    HighFive::DataSet dmat = file.createDataSet<double>(
          "/4elCouplingMat", HighFive::DataSpace(dims));

    // start point and the number of the total number of elements
    // to be written by this process
    size_t start = mpi->divideWorkIter(numTriplets)[0] * bandProd;
    size_t offset = start;

    // Note: HDF5 < v1.10.2 cannot write datasets larger than 2 Gbs
    // ( due to max(int 32 bit))/1024^3 = 2Gb overflowing in MPI)
    // In order to be compatible with older versions, we split the tensor
    // into smaller chunks and write them to separate datasets
    // slower, but it will work more often.

    // maxSize represents 2GB worth of std::complex<doubles>, since that's
    // what we write
    auto maxSize = int(pow(1000, 3)) / sizeof(double);
    size_t smallestSize = bandProd; // 1 point pair
    std::vector<int> bunchSizes;

    // determine the size of each bunch of electronic bravais vectors
    // the BunchSizes vector tells us how many are in each set
    int numTripletsBunch = mpi->divideWorkIter(numTriplets).back() + 1 - mpi->divideWorkIter(numTriplets)[0];

    int bunchSize = 0;
    for (int i = 0; i < numTripletsBunch; i++) {
      bunchSize++;
      // this bunch is as big as possible, stop adding to it
      if ((bunchSize + 1) * smallestSize > maxSize) {
        bunchSizes.push_back(bunchSize);
        bunchSize = 0;
      }
    }
    // push the last one, no matter the size, to the list of bunch sizes
    bunchSizes.push_back(bunchSize);

    // determine the number of bunches. The last bunch may be smaller than the rest
    int numDatasets = bunchSizes.size();

    // we now loop over these data sets and write each chunk in parallel
    int netOffset = 0;  // offset from first bunch in this set to current bunch

    for (int iBunch = 0; iBunch < numDatasets; iBunch++) {

      // we need to determine the start, stop and offset of this
      // sub-slice of the dataset available to this process
      size_t bunchElements = bunchSizes[iBunch] * smallestSize;
      size_t bunchOffset = offset + netOffset;
      netOffset += bunchElements;

      // Each process writes to hdf5
      // The format is ((startRow,startCol),(numRows,numCols)).write(data)
      // Because it's a vector (1 row) all processes write to row=0, col=startPoint
      // with nRows = 1, nCols = number of items this process will write.
      dmat.select({0, bunchOffset}, {1, bunchElements}).write_raw(&allMatrixElementsSq[0]);

    }
    } // end parallel write section
    #else
    {

    // throw an error if there are too many elements to write
    unsigned int globalSize = numTriplets * bandProd;
    auto maxSize = int(pow(1000, 3)) / sizeof(double);
    if(globalSize > maxSize) {
      Error("Your requested el-el matrix element file size is greater than the allowed size\n"
        "for a single write by HDF5. Either compile Phoebe with a parallel copy of HDF5 or\n"
        "request to output fewer matrix elements. ");
    }

    // call an mpi collective to gather matrix elements
    std::vector<double> collectedGs;
    if(mpi->mpiHead()) collectedGs.resize(globalSize);
    mpi->allGatherv(&collectedGs, &allMatrixElementsSq);

    // write elph matrix elements
    HighFive::File file(outFileName, HighFive::File::Overwrite);
    file.createDataSet("/4elCouplingMat", collectedGs);

    }
    #endif
    // now we write a few other pieces of smaller information using only mpiHead
    if (mpi->mpiHead()) {

      HighFive::File file(outFileName, HighFive::File::ReadWrite);

      // shape the points pairs list into a format that can be written
      Eigen::MatrixXd pointsTemp(pointsTriplets.size(),9);
      for (size_t iTriplet = 0; iTriplet < pointsTriplets.size(); iTriplet++) {

        auto thisTriplet = pointsTriplets[iTriplet];
        Eigen::Vector3d k1C = kPoints.cartesianToCrystal(std::get<0>(thisTriplet));
        Eigen::Vector3d k2C = kPoints.cartesianToCrystal(std::get<1>(thisTriplet));
        Eigen::Vector3d k3C = kPoints.cartesianToCrystal(std::get<2>(thisTriplet));
        for( int i : {0,1,2} ) {
          pointsTemp(iTriplet,i) = k1C(i);
          pointsTemp(iTriplet,i+3) = k2C(i);
          pointsTemp(iTriplet,i+6) = k3C(i);
        }
      }
      // write the points pairs to file
      file.createDataSet("/pointsTriplets", pointsTemp);

      // write the band ranges
      std::vector<int> temp;
      temp.push_back(g2PlotEl1Bands.first);
      temp.push_back(g2PlotEl1Bands.second);
      file.createDataSet("/elBandRange1", temp);
      temp.clear();
      temp.push_back(g2PlotEl2Bands.first);
      temp.push_back(g2PlotEl2Bands.second);
      file.createDataSet("/elBandRange2", temp);
      temp.clear();
      temp.push_back(g2PlotEl3Bands.first);
      temp.push_back(g2PlotEl3Bands.second);
      file.createDataSet("/elBandRange3", temp);
      temp.clear();
      temp.push_back(g2PlotEl4Bands.first);
      temp.push_back(g2PlotEl4Bands.second);
      file.createDataSet("/elBandRange4", temp);

    }
  } catch (std::exception &error) {
      Error("Issue writing el-el Wannier representation to hdf5.");
  }
  // close else for HDF5_AVAIL
  #else
  Error("You cannot output the elph matrix elements to HDF5 because your copy of \n"
        "Phoebe has not been compiled with HDF5 support.");
  #endif
}

//TODO
// * check that the bands supplied by users don't exceed nWannier

void ElElCouplingPlotApp::checkRequirements(Context &context) {
  throwErrorIfUnset(context.getElectronH0Name(), "electronH0Name");
  //if(context.getG2PlotStyle() == "pointsPath")
  //  throwErrorIfUnset(context.getPathExtrema(), "points path extrema");
  throwErrorIfUnset(context.getKMesh(), "kMesh");
  // TODO this does not work!
  throwErrorIfUnset(context.getG2PlotFixedPoint(), "couplingFixedPoint");
  // check that crystal structure was provided
  std::string crystalMsg = "crystal structure";
  throwErrorIfUnset(context.getInputAtomicPositions(), crystalMsg);
  throwErrorIfUnset(context.getInputSpeciesNames(), crystalMsg);
  throwErrorIfUnset(context.getInputAtomicSpecies(), crystalMsg);
}
