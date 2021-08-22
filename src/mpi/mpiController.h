#ifndef MPI_CONTROLLER_H
#define MPI_CONTROLLER_H

#include <chrono>
#include <complex>
#include <vector>
#include "eigen.h"
#include <tuple>
#include "exceptions.h"
#include <Kokkos_Core.hpp>

#ifdef MPI_AVAIL
#include <mpi.h>
#endif

const int worldComm_ = 0;
const int intraPoolComm_ = 1;
const int interPoolComm_ = 2;

/* NOTE: When using this object make sure to use the divideWork
functions to set up the initial division of tasks --
functions here (see gather) operate under the
assumption that they will gather using the same distribution of
 labor as divideWork provides.*/

/** Class for handling the MPI library usage inside of phoebe.
 * We define 3 communicators.
 * 1) MPI_COMM_WORLD: this is the communicator involving all MPI processes
 *
 * If the -ps flag is used on the command line, the grid of MPI processes is
 * viewed as a rectangle, and we define two MPI groups:
 * 2) intraPoolComm: the rows of the rectangle of MPI processes identify pools.
 * The number of processes (columns) in each pool is the number specified
 * with the `ps` flag from the command line.
 * 3) interPoolComm: is a communicator grouping the columns of the MPI rectangle
 * built after specifying the `ps` flag.
 * Note that in order for the `ps` flag to work, we require that all MPI
 * processes in the world communicator can be divided in the rectangle with no
 * MPI process left idle.
 */
class MPIcontroller {
 private:

  // default comm "MPI_COMM_WORLD"
  int size = 0;  // number of MPI processes
  int rank;
  const int mpiHeadId = 0;
  const int mpiHeadPoolId = 0;
  const int mpiHeadColsId = 0;

  int poolSize = 1; // # of MPI processes in the pool
  int poolRank = 0; // rank of the MPI process within the pool from 0 to poolSize
  int poolId = 0; // id of the pool
#ifdef MPI_AVAIL
  MPI_Comm intraPoolCommunicator;
  MPI_Comm interPoolCommunicator;
  MPI_Comm worldCommunicator = MPI_COMM_WORLD;
#endif

  // helper function used internally
  std::tuple<std::vector<int>, std::vector<int>> workDivHelper(size_t numTasks) const;

#ifdef MPI_AVAIL
  double startTime;  // the time for the entire mpi operation
  /** Utility function used to convert the integer communicator to a MPI_COMM
   * communicator.
   * @param communicator: a MPI_COMM value of th communicator.
   * @return
   */
  std::tuple<MPI_Comm, int> decideCommunicator(const int& communicator) const;
#else
  std::chrono::steady_clock::time_point startTime;
#endif

 public:
  // MPIcontroller class constructors -----------------------------------
  /** a constructor which sets up the MPI environment, initializes the
   * communicator, and starts a timer **/
  MPIcontroller(int argc, char *argv[]);

  /** Calls finalize and potentially reports statistics */
  void finalize() const;

  // Collective communications functions -----------------------------------
  /** Wrapper for the MPI_Broadcast function.
   *  @param dataIn: pointer to data structure to broadcast
   */
  template <typename T>
  void bcast(T* dataIn, const int& communicator=worldComm) const;

  /** Wrapper for MPI_Reduce in the case of a summation.
   * @param dataIn: pointer to sent data from each rank.
   * @param dataOut: pointer to buffer to receive summed data.
   */
  template <typename T>
  void allReduceSum(T* dataIn, T* dataOut) const;

  /** Wrapper for MPI_AllReduce in the case of a summation in-place.
   * @param data: pointer to sent data from each rank.
   * Gets overwritten with the result of the MPI allreduce('SUM') operation.
   */
  template <typename T>
  void allReduceSum(T* dataIn, const int& communicator=worldComm) const;

  /** Wrapper for MPI_Reduce in the case of a summation.
   * @param dataIn: pointer to sent data from each rank,
   *       also acts as a receive buffer, as reduce is implemented IP.
   */
  template <typename T>
  void reduceSum(T* dataIn) const;

  /** Wrapper for MPI_Reduce which identifies the maximum of distributed data
   * @param dataIn: pointer to sent data from each rank.
   *       also acts as a receive buffer, as reduce is implemented IP.
   */
  template <typename T>
  void reduceMax(T* dataIn) const;

  /** Wrapper for MPI_AllReduce which identifies the maximum of distributed data
   * @param dataIn: pointer to sent data from each rank.
   *       also acts as a receive buffer, as reduce is implemented IP.
   */
  template <typename T>
  void allReduceMax(T* dataIn, const int& communicator=worldComm) const;

  /** Wrapper for MPI_Reduce which identifies the minimum of distributed data
   * @param dataIn: pointer to sent data from each rank.
   *       also acts as a receive buffer, as reduce is implemented IP.
   */
  template <typename T>
  void reduceMin(T* dataIn) const;

  /** Wrapper for MPI_AllReduce which identifies the minimum of distributed data
   * @param dataIn: pointer to sent data from each rank.
   *       also acts as a receive buffer, as reduce is implemented IP.
   */
  template <typename T>
  void allReduceMin(T* dataIn) const;

  /** Wrapper for MPI_Gatherv which collects data from different ranks
   * (with the possibility of a different number of elements from each
   * process) and combines it into one buffer.
   * @param dataIn: pointer to sent data from each rank, with length
   *       of the number of points belonging to this rank.
   * @param dataOut: pointer to output buffer, allocated only by the
   *       head rank, of length to contain all data from all processes.
   */
  template <typename T, typename V>
  void gatherv(T* dataIn, V* dataOut) const;

  /** Wrapper for MPI_Gather which collects data from different ranks
   * and combines it into one buffer.
   * @param dataIn: pointer to sent data from each rank, with length
   *       of the number of points belonging to this rank.
   * @param dataOut: pointer to output buffer, allocated only by the
   *       head rank, of length to contain all data from all processes.
   */
  template <typename T, typename V>
  void gather(T* dataIn, V* dataOut) const;

  /** Wrapper for MPI_Allgatherv which collects data from different ranks
   * (with the possibility of a different number of elements from each
   * process) and combines it into one buffer, which is also broadcast
   * to all processes.
   * @param dataIn: pointer to sent data from each rank, with length
   *       of the number of points belonging to this rank.
   * @param dataOut: pointer to output buffer, allocated only by the
   *       head rank, of length to contain all data from all processes.
   */
  template <typename T, typename V>
  void allGatherv(T* dataIn, V* dataOut) const;

  /** Wrapper for MPI_Allgather which collects data from different ranks
   * and combines it into one buffer, which is also broadcast
   * to all processes.
   * @param dataIn: pointer to sent data from each rank, with length
   *       of the number of points belonging to this rank.
   * @param dataOut: pointer to output buffer, allocated only by the
   *       head rank, of length to contain all data from all processes.
   */
  template <typename T, typename V>
  void allGather(T* dataIn, V* dataOut) const;

  template <typename T>
  void datatypeHelper(MPI_Datatype* container, MPI_Count count, T* data) const;

  template <typename T>
  void bigAllGatherV(T* dataIn, T* dataOut, std::vector<size_t> workDivs_, std::vector<size_t> workDivisionHeads_) const;

  // Asynchronous functions
  /** Wrapper for MPI_Barrier()
   */
  void barrier() const;

  // Utility functions -----------------------------------
  /** Simple function to tell us if this process is the head
   * @return isRank: returns true if this rank is the head.
   */
  bool mpiHead() const { return rank == mpiHeadId; }

  /** Returns true if the MPI process belongs to the head pool
   * @return isHeadPool: returns true if this MPI process is in the head pool.
   */
  bool mpiHeadPool() const { return poolId == mpiHeadPoolId; }

  /** Function to return the rank of a process.
   * @return rank: the rank of this process.
   */
  int getRank(const int& communicator=worldComm) const {
    if (communicator == worldComm) {
      return rank;
    } else if (communicator == intraPoolComm) {
      return poolRank;
    } else {
      Error("Invalid communicator in getSize");
      return 0;
    }
  };

  /** Function to return the number of ranks available.
   * @return size: number of ranks
   */
  int getSize(const int& communicator=worldComm) const {
    if (communicator == worldComm) {
      return size;
    } else if (communicator == intraPoolComm) {
      return poolSize;
    } else {
      Error("Invalid communicator in getSize");
      return 0;
    }
  };

  // Error reporting and statistics
  void errorReport(int errCode) const;  // collect errors from processes and
                                        // reports them, then kills the code

  /** Function to return the time elapsed since the mpi environment started.
  */
  void time() const;

  /** Divides a number of tasks appropriately for the current MPI env.
   * @return divs: returns a vector of length 2, containing start and stop
   *       points for the divided number of tasks.
   */
  std::vector<size_t> divideWork(size_t numTasks);  // divide up a set of work

  /** Divides a number of tasks appropriately for the current MPI env.
   * @return divs: returns an iterator of points for the divided number of tasks.
   */
  std::vector<size_t> divideWorkIter(size_t numTasks, const int& communicator=worldComm);

  /** integer used to specify the call to MPI uses the world communicator.
   */
  static const int worldComm;

  /** integer used to specify the call to MPI uses the pool communicator.
   */
  static const int intraPoolComm;

  /** integer used to specify the call to MPI uses the inter-Pool communicator.
   */
  static const int interPoolComm;
};

// we need to use the concept of a "type traits" object to serialize the
// standard cpp types and then we can define specific implementations of this
// for types which are not standard
namespace mpiContainer {

#ifdef MPI_AVAIL
// Forward declaration for a basic container type
template<typename...>
struct containerType;

// Define a macro to shorthand define this container for scalar types
// Size for basic scalar types will always be 1
#define MPIDataType(cppType, mpiType)                                 \
  template<>                                                          \
  struct containerType<cppType> {                                     \
    static inline cppType *getAddress(cppType *data) { return data; } \
    static inline size_t getSize(cppType *data) {                     \
      if (!data) return 0;                                            \
      return 1;                                                       \
    }                                                                 \
    static inline MPI_Datatype getMPItype() { return mpiType; }       \
  };

  // Use definition to generate containers for scalar types
  MPIDataType(int, MPI_INT)
  MPIDataType(long, MPI_LONG)
  MPIDataType(unsigned int, MPI_UNSIGNED)
  MPIDataType(float, MPI_FLOAT)
  MPIDataType(double, MPI_DOUBLE)
  MPIDataType(std::complex<double>, MPI_DOUBLE_COMPLEX)
  MPIDataType(std::complex<float>, MPI_COMPLEX)
  MPIDataType(size_t, MPI_UNSIGNED_LONG_LONG)

#undef MPIDataType

    // A container for a std::vector
  template<typename T>
  struct containerType<std::vector<T>> {
    static inline T *getAddress(std::vector<T> *data) { return data->data(); }
    static inline size_t getSize(std::vector<T> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for <Eigen::Matrix<T, -1, 1, 0, -1, 1>
  template<typename T>
  struct containerType<Eigen::Matrix<T, -1, 1, 0, -1, 1>> {
    static inline T *getAddress(Eigen::Matrix<T, -1, 1, 0, -1, 1> *data) { return data->data(); }
    static inline size_t getSize(Eigen::Matrix<T, -1, 1, 0, -1, 1> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for <Eigen::Matrix<T, -1, -1>
  template<typename T>
  struct containerType<Eigen::Matrix<T, -1, -1>> {
    static inline T *getAddress(Eigen::Matrix<T, -1, -1> *data) { return data->data(); }
    static inline size_t getSize(Eigen::Matrix<T, -1, -1> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for Eigen::Tensor<T, 5>
  template<typename T>
  struct containerType<Eigen::Tensor<T, 5>> {
    static inline T *getAddress(Eigen::Tensor<T, 5> *data) { return data->data(); }
    static inline size_t getSize(Eigen::Tensor<T, 5> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for Eigen::Tensor<T, 4>
  template<typename T>
  struct containerType<Eigen::Tensor<T, 4>> {
    static inline T *getAddress(Eigen::Tensor<T, 4> *data) { return data->data(); }
    static inline size_t getSize(Eigen::Tensor<T, 4> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for Eigen::Tensor<T, 3>
  template<typename T>
  struct containerType<Eigen::Tensor<T, 3>> {
    static inline T *getAddress(Eigen::Tensor<T, 3> *data) { return data->data(); }
    static inline size_t getSize(Eigen::Tensor<T, 3> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype(); }
  };
  // Container for Eigen::MatrixXi
  template<>
  struct containerType<Eigen::MatrixXi> {
    static inline int *getAddress(Eigen::MatrixXi *data) { return data->data(); }
    static inline size_t getSize(Eigen::MatrixXi *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype(); }
  };
  // Container for Eigen::VectorXi
  template<>
  struct containerType<Eigen::VectorXi> {
    static inline int *getAddress(Eigen::VectorXi *data) { return data->data(); }
    static inline size_t getSize(Eigen::VectorXi *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype(); }
  };
  // Container for Eigen::Vector3i
  template<>
  struct containerType<Eigen::Vector3i> {
    static inline int *getAddress(Eigen::Vector3i *data) { return data->data(); }
    static inline size_t getSize(Eigen::Vector3i *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype(); }
  };
  // Container for Eigen::Vector3d
  template<>
  struct containerType<Eigen::Vector3d> {
    static inline double *getAddress(Eigen::Vector3d *data) { return data->data(); }
    static inline size_t getSize(Eigen::Vector3d *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<double>::getMPItype(); }
  };
  // Container for Eigen::VectorXcd
  template<>
  struct containerType<Eigen::VectorXcd> {
    static inline std::complex<double> *getAddress(Eigen::VectorXcd *data) { return data->data(); }
    static inline size_t getSize(Eigen::VectorXcd *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<std::complex<double>>::getMPItype(); }
  };
  // Container for Eigen::VectorXd
  template<>
  struct containerType<Eigen::VectorXd> {
    static inline double *getAddress(Eigen::VectorXd *data) { return data->data(); }
    static inline size_t getSize(Eigen::VectorXd *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return containerType<double>::getMPItype(); }
  };

  template<>
  struct containerType<Kokkos::View<Kokkos::complex<double>****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>> {
    static inline Kokkos::complex<double> *getAddress(Kokkos::View<Kokkos::complex<double>****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> *data) { return data->data(); }
    static inline size_t getSize(Kokkos::View<Kokkos::complex<double>****, Kokkos::LayoutRight, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>> *data) { return data->size(); }
    static inline MPI_Datatype getMPItype() { return MPI_COMPLEX16; }
  };
#endif
}  // namespace mpiContainer

// Collective communications functions -----------------------------------
template <typename T>
void MPIcontroller::bcast(T* dataIn, const int& communicator) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  if (communicator == intraPoolComm && poolSize == 1) return;

  auto t = decideCommunicator(communicator);
  MPI_Comm comm = std::get<0>(t);
  int broadcasterId = std::get<1>(t);

  int errCode = MPI_Bcast(containerType<T>::getAddress(dataIn),
                      containerType<T>::getSize(dataIn),
                      containerType<T>::getMPItype(), broadcasterId, comm);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::reduceSum(T* dataIn) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;

  if (rank == 0) {
    errCode =
        MPI_Reduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                   containerType<T>::getSize(dataIn),
                   containerType<T>::getMPItype(), MPI_SUM, mpiHeadId, MPI_COMM_WORLD);
  } else {
    errCode = MPI_Reduce(
        containerType<T>::getAddress(dataIn),
        containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
        containerType<T>::getMPItype(), MPI_SUM, mpiHeadId, MPI_COMM_WORLD);
  }
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::allReduceSum(T* dataIn, T* dataOut) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;

  errCode = MPI_Allreduce(
      containerType<T>::getAddress(dataIn),
      containerType<T>::getAddress(dataOut), containerType<T>::getSize(dataIn),
      containerType<T>::getMPItype(), MPI_SUM, MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::reduceMax(T* dataIn) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;

  if (rank == 0) {
    errCode =
        MPI_Reduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                   containerType<T>::getSize(dataIn),
                   containerType<T>::getMPItype(), MPI_MAX, mpiHeadId, MPI_COMM_WORLD);
  } else {
    errCode = MPI_Reduce(
        containerType<T>::getAddress(dataIn),
        containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
        containerType<T>::getMPItype(), MPI_MAX, mpiHeadId, MPI_COMM_WORLD);
  }
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::allReduceSum(T* dataIn, const int& communicator) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  if (communicator == intraPoolComm && poolSize == 1) return;

  MPI_Comm comm = std::get<0>(decideCommunicator(communicator));

  int errCode =
      MPI_Allreduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                    containerType<T>::getSize(dataIn),
                    containerType<T>::getMPItype(), MPI_SUM, comm);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::allReduceMax(T* dataIn, const int& communicator) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  if (size == 1) return;
  if (communicator == intraPoolComm && poolSize == 1) return;

  MPI_Comm comm = std::get<0>(decideCommunicator(communicator));

  int errCode =
      MPI_Allreduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                    containerType<T>::getSize(dataIn),
                    containerType<T>::getMPItype(), MPI_MAX, comm);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #endif
}

template <typename T>
void MPIcontroller::reduceMin(T* dataIn) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;

  if (rank == 0) {
    errCode =
        MPI_Reduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                   containerType<T>::getSize(dataIn),
                   containerType<T>::getMPItype(), MPI_MIN, mpiHeadId, MPI_COMM_WORLD);
  } else {
    errCode = MPI_Reduce(
        containerType<T>::getAddress(dataIn),
        containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
        containerType<T>::getMPItype(), MPI_MIN, mpiHeadId, MPI_COMM_WORLD);
  }
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

template <typename T>
void MPIcontroller::allReduceMin(T* dataIn) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;
  errCode =
      MPI_Allreduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                    containerType<T>::getSize(dataIn),
                    containerType<T>::getMPItype(), MPI_MIN, MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #endif
}

/* ---------- gather function wrappers ------------- */

template <typename T>
void pointerSwap(T* dataIn, std::vector<T>* dataOut) {
  std::fill(dataOut->begin(), dataOut->end(), *dataIn);
}

template <typename T> // this is implemented only for std::vector
void pointerSwap(T* dataIn, T* dataOut) {
  for (unsigned int i=0;i<dataIn->size();i++) {
    (*dataOut)[i] = (*dataIn)[i];
  }
}

template <typename T, typename V>
void MPIcontroller::gatherv(T* dataIn, V* dataOut) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  int errCode;

  // calculate the number of elements coming from each process
  // this will correspond to the save division of elements
  // as divideWorkIter provides.
  int numTasks = containerType<V>::getSize(dataOut);
  auto tup = workDivHelper(numTasks);
  std::vector<int> workDivs = std::get<0>(tup);
  std::vector<int> workDivisionHeads = std::get<1>(tup);

  errCode = MPI_Gatherv(
      containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
      containerType<T>::getMPItype(), containerType<V>::getAddress(dataOut),
      workDivs.data(), workDivisionHeads.data(), containerType<V>::getMPItype(),
      mpiHeadId, MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #else
  pointerSwap(dataIn, dataOut);  // just switch the pointers in serial case
  #endif
}

template <typename T, typename V>
void MPIcontroller::gather(T* dataIn, V* dataOut) const {
    using namespace mpiContainer;
  #ifdef MPI_AVAIL
  int errCode;

  errCode = MPI_Gather(
      containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
      containerType<T>::getMPItype(), containerType<V>::getAddress(dataOut),
      containerType<T>::getSize(dataIn), containerType<V>::getMPItype(),
      mpiHeadId, MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #else
  pointerSwap(dataIn, dataOut);  // just switch the pointers in serial case
  #endif
}


template <typename T, typename V>
void MPIcontroller::allGatherv(T* dataIn, V* dataOut) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  int errCode;

  // calculate the number of elements coming from each process
  // this will correspond to the save division of elements
  // as divideWorkIter provides.
  size_t numTasks = containerType<V>::getSize(dataOut);
  auto tup = workDivHelper(numTasks);
  std::vector<int> workDivs = std::get<0>(tup);
  std::vector<int> workDivisionHeads = std::get<1>(tup);

  errCode = MPI_Allgatherv(
      containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
      containerType<T>::getMPItype(), containerType<V>::getAddress(dataOut),
      workDivs.data(), workDivisionHeads.data(), containerType<V>::getMPItype(),
      MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #else
  pointerSwap(dataIn, dataOut);
  #endif
}

template <typename T, typename V>
void MPIcontroller::allGather(T* dataIn, V* dataOut) const {
    using namespace mpiContainer;
  #ifdef MPI_AVAIL
  int errCode;

  errCode = MPI_Allgather(
      containerType<T>::getAddress(dataIn), containerType<T>::getSize(dataIn),
      containerType<T>::getMPItype(), containerType<V>::getAddress(dataOut),
      containerType<T>::getSize(dataIn), containerType<V>::getMPItype(),
      MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
  #else
  pointerSwap(dataIn, dataOut);  // just switch the pointers in serial case
  #endif
}

template <typename T>
void MPIcontroller::datatypeHelper(MPI_Datatype* container, MPI_Count count, T* data) const {

 using namespace mpiContainer;
 #ifdef MPI_AVAIL

    // these hold the count of elements sent in each block, and
    // the remainder in case the division is not even among the ranks
    MPI_Count blockSize = count / INT_MAX;
    MPI_Count remain = count % INT_MAX;

    std::cout << "blockSize remain count " << blockSize << " " << remain << " " << count << std::endl;

    /* first, we create two intermediate data types, block and remainder.
    *
    *  MPI_TYPE_VECTOR: replication of a datatype into locations that consist of equally spaced blocks
    *  MPI_TYPE_VECTOR(COUNT, BLOCKLENGTH, STRIDE, OLDTYPE, NEWTYPE, IERROR)
    *  blocks of data will be divided into INT_MAX chunks */
    MPI_Datatype block;
    MPI_Type_vector(blockSize, INT_MAX, INT_MAX, containerType<T>::getMPItype(), &block);

    // set up a container for the remainder of the data as well
    MPI_Datatype remainder;
    MPI_Type_contiguous(remain, containerType<T>::getMPItype(), &remainder);

    // need to set up the "extent" of the data
    MPI_Aint lb, extent;
    // This function gets the lower bound and extent for a datatype
    MPI_Type_get_extent(containerType<T>::getMPItype(), &lb, &extent);

    // Find address displacements of data
    MPI_Aint remdisp          = (MPI_Aint)blockSize * INT_MAX * extent;
    MPI_Aint offset       = 0;
    int blocklengths[2]       = {1,1};
    MPI_Aint displacements[2] = {offset, offset + remdisp};

    // create the datatype to actually be used,
    // constructed from these block and remainder types
    MPI_Datatype types[2]     = {block, remainder};
    MPI_Type_create_struct(2, blocklengths, displacements, types, container);

    std::cout << "just before freeing block and remainder" << std::endl;

    // free the intermediate datatypes
    MPI_Type_free(&block);
    MPI_Type_free(&remainder);

    std::cout << "just after freeing block and remainder " << std::endl;
    // establish the container type within MPI
    MPI_Type_commit(container);

    std::cout << "commited container "<< std::endl;

  #endif
}


template <typename T>
void MPIcontroller::bigAllGatherV(T* dataIn, T* dataOut, std::vector<size_t> workDivs_, std::vector<size_t> workDivisionHeads_) const {

  using namespace mpiContainer;
  #ifdef MPI_AVAIL

    // if there's only one process we don't need to act
   // if (size == 1) return;
    int errCodeSend, errCodeRecv;

    // calculate the number of elements coming from each process
    // this will correspond to the save division as divideWorkIter
    std::vector<MPI_Aint> workDivisionHeads(size);
    std::vector<MPI_Count> workDivs(size);

    workDivisionHeads[0] = (MPI_Aint)workDivisionHeads_[0];
    workDivs[0] = (MPI_Count)workDivs_[0];
    if(rank == 0) std::cout << "rank, workDivs#, workDivisionHeads# " << rank << " " << workDivs[0] << " " << workDivisionHeads[0] << std::endl;

    for (int i = 1; i < size; i++) {
      workDivisionHeads[i] = (MPI_Aint)(workDivisionHeads_[i]);
      workDivs[i] = (MPI_Count)workDivs_[i];
      if(rank == 0) std::cout << "rank, workDivs#, workDivisionHeads# " << i << " " << workDivs[i] << " " << workDivisionHeads[i] << std::endl;
    }
    MPI_Aint * recvdispls = workDivisionHeads.data();

    // this is required to use non-blocking communications,
    // which will help this method be a bit faster
    std::vector<MPI_Request> reqs(2*size);

    // need this to determine where the received chunks will be located
    // in the final dataOut buffer
    MPI_Aint lb, recvextent;
    MPI_Type_get_extent(containerType<T>::getMPItype(), &lb, &recvextent);

    // this process receives from all other processes
    for (int i=0; i<size; i++) {

      MPI_Datatype container;
      datatypeHelper(&container, workDivs[i], dataOut);

      std::cout << "rank i recvCount dataOut.address reqs.size " << rank << " "<< i << " " << workDivs[i] << " " << containerType<T>::getAddress(dataOut) << " " << reqs.size() << std::endl;

      // MPI_Irecv is a non-blocking communication method
      errCodeRecv = MPI_Irecv(containerType<T>::getAddress(dataOut)+recvdispls[i]*recvextent,
          1, container, i, 0, MPI_COMM_WORLD, &reqs[i]);

      MPI_Type_free(&container);
      if (errCodeRecv != MPI_SUCCESS) errorReport(errCodeRecv);

    }

    // send the data from this process to all other processes
    for (int j=rank; j<(size+rank); j++) {

      // smart way to balance communications
      int i = j%size;

      MPI_Datatype container;
      datatypeHelper(&container, workDivs[rank], dataIn);

      //if(rank != 0)
      std::cout << "rank i sendCount dataIn.address reqs.size size+i " << rank << " "<< i << " " << workDivs[rank] << " " << containerType<T>::getAddress(dataIn) << " " << reqs.size() << " " << size+i << std::endl;

      errCodeSend = MPI_Isend(containerType<T>::getAddress(dataIn), 1, container,
          i, 0, MPI_COMM_WORLD, &reqs[size+i]);
      std::cout << " send from rank " << rank << " got through." << std::endl;

      MPI_Type_free(&container);
      if (errCodeSend != MPI_SUCCESS) errorReport(errCodeSend);
    }

    MPI_Waitall(2*size, reqs.data(), MPI_STATUSES_IGNORE);

  #endif
}

#endif
