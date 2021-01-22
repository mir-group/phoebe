#ifndef MPICONTROLLER_H
#define MPICONTROLLER_H

#include <chrono>
#include <complex>
#include <vector>
#include "eigen.h"
#include <tuple>
#include "exceptions.h"

#ifdef MPI_AVAIL
#include <mpi.h>
#endif

/* NOTE: When using this object make sure to use the divideWork
functions to set up the initial division of tasks --
functions here (see gather) operate under the
assumption that they will gather using the same distribution of
 labor as divideWork provides.*/

/** Class for handling the MPI library usage inside of phoebe.
 */
class MPIcontroller {
 private:

  // default comm "MPI_COMM_WORLD"
  int size = 0;  // number of MPI processses
  int rank;
  const int mpiHeadId = 0;

  // helper function used internally
  std::tuple<std::vector<int>,std::vector<int>> workDivHelper(size_t numTasks) const;

  #ifdef MPI_AVAIL
    double startTime;  // the time for the entire mpi operation
  #else
    std::chrono::steady_clock::time_point startTime;
  #endif

 public:
  // MPIcontroller class constructors -----------------------------------
  /** a constructor which sets up the MPI environment, initializes the
   * communicator, and starts a timer **/
  MPIcontroller();

  /** Calls finalize and potentially reports statistics */
  void finalize() const;

  // Collective communications functions -----------------------------------
  /** Wrapper for the MPI_Broadcast function.
   *  @param dataIn: pointer to data structure to broadcast
   */
  template <typename T>
  void bcast(T* dataIn) const;

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
  void allReduceSum(T* dataIn) const;

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
  void allReduceMax(T* dataIn) const;

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

  // point to point functions -----------------------------------
  // template<typename T> void send(T&& data) const;
  // template<typename T> void recv(T&& data) const;

  // Asynchronous functions
  /** Wrapper for MPI_Barrier()
   */
  void barrier() const;

  // Utility functions -----------------------------------
  /** Simple function to tell us if this process is the head
   * @return isRank: returns true if this rank is the head.
   */
  bool mpiHead() const { return rank == mpiHeadId; }
  /** Function to return the rank of a process.
   * @return rank: the rank of this process.
   */
  int getRank() const { return rank; }
  /** Function to return the number of ranks available.
   * @return size: number of ranks
   */
  int getSize() const { return size; }

  // Error reporting and statistics
  void errorReport(int errCode) const;  // collect errors from processes and
                                        // reports them, then kills the code
  void time() const;  // returns time elapsed since mpi started

  /** Divides a number of tasks appropriately for the current MPI env.
   * @return divs: returns a vector of length 2, containing start and stop
   *       points for the divided number of tasks.
   */
  std::vector<int> divideWork(size_t numTasks);  // divide up a set of work
  /** Divides a number of tasks appropriately for the current MPI env.
   * @return divs: returns an iterator of points for the divided number of tasks.
   */
  std::vector<int> divideWorkIter(size_t numTasks);

};

// we need to use the concept of a "type traits" object to serialize the
// standard cpp types and then we can define specific implementations of this
// for types which are not standard
namespace mpiContainer {
        #ifdef MPI_AVAIL
        // Forward declaration for a basic container type
        template <typename...>
        struct containerType;

        // Define a macro to shorthand define this container for scalar types
        // Size for basic scalar types will always be 1
        #define MPIDataType(cppType, mpiType)                                 \
          template <>                                                         \
          struct containerType<cppType> {                                     \
            static inline cppType* getAddress(cppType* data) { return data; } \
            static inline size_t getSize(cppType* data) {                     \
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

        #undef MPIDataType

        // A container for a std::vector
        template <typename T> struct containerType<std::vector<T>> {
                static inline T* getAddress(std::vector<T>* data) { return data->data(); }
                static inline size_t getSize(std::vector<T>* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype();}
        };
        // Container for <Eigen::Matrix<T, -1, 1, 0, -1, 1>
        template <typename T> struct containerType<Eigen::Matrix<T, -1, 1, 0, -1, 1>> {
                static inline T* getAddress(Eigen::Matrix<T, -1, 1, 0, -1, 1>* data) { return data->data(); }
                static inline size_t getSize(Eigen::Matrix<T, -1, 1, 0, -1, 1>* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype();}
        };
        // Container for <Eigen::Matrix<T, -1, -1>
        template <typename T> struct containerType<Eigen::Matrix<T, -1, -1>> {
                static inline T* getAddress(Eigen::Matrix<T, -1, -1>* data) { return data->data(); }
                static inline size_t getSize(Eigen::Matrix<T, -1, -1>* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype();}
        };
        // Container for Eigen::Tensor<T, 5>
        template <typename T> struct containerType<Eigen::Tensor<T, 5>> {
                static inline T* getAddress(Eigen::Tensor<T, 5>* data) { return data->data(); }
                static inline size_t getSize(Eigen::Tensor<T, 5>* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype();}
        };
        // Container for Eigen::Tensor<T, 3>
        template <typename T> struct containerType<Eigen::Tensor<T, 3>> {
                static inline T* getAddress(Eigen::Tensor<T, 3>* data) { return data->data(); }
                static inline size_t getSize(Eigen::Tensor<T, 3>* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<T>::getMPItype();}
        };
        // Container for Eigen::MatrixXi
        template <> struct containerType<Eigen::MatrixXi> {
                static inline int* getAddress(Eigen::MatrixXi* data) { return data->data(); }
                static inline size_t getSize(Eigen::MatrixXi* data) { return data->size(); }
                static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype();}
        };
        // Container for Eigen::VectorXi
        template <> struct containerType<Eigen::VectorXi> {
          static inline int* getAddress(Eigen::VectorXi* data) { return data->data(); }
          static inline size_t getSize(Eigen::VectorXi* data) { return data->size(); }
          static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype();}
        };
        // Container for Eigen::Vector3i
        template <> struct containerType<Eigen::Vector3i> {
          static inline int* getAddress(Eigen::Vector3i* data) { return data->data(); }
          static inline size_t getSize(Eigen::Vector3i* data) { return data->size(); }
          static inline MPI_Datatype getMPItype() { return containerType<int>::getMPItype();}
        };
        // Container for Eigen::VectorXcd
        template <> struct containerType<Eigen::VectorXcd> {
          static inline std::complex<double>* getAddress(Eigen::VectorXcd* data) { return data->data(); }
          static inline size_t getSize(Eigen::VectorXcd* data) { return data->size(); }
          static inline MPI_Datatype getMPItype() { return containerType<std::complex<double>>::getMPItype();}
        };

#endif
}  // namespace mpiContainer

// Collective communications functions -----------------------------------
template <typename T>
void MPIcontroller::bcast(T* dataIn) const {
  using namespace mpiContainer;
#ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;

  errCode = MPI_Bcast(containerType<T>::getAddress(dataIn),
                      containerType<T>::getSize(dataIn),
                      containerType<T>::getMPItype(), mpiHeadId, MPI_COMM_WORLD);
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
void MPIcontroller::allReduceSum(T* dataIn) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;
  errCode =
      MPI_Allreduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                    containerType<T>::getSize(dataIn),
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
void MPIcontroller::allReduceMax(T* dataIn) const {
  using namespace mpiContainer;
  #ifdef MPI_AVAIL
  if (size == 1) return;
  int errCode;
  errCode =
      MPI_Allreduce(MPI_IN_PLACE, containerType<T>::getAddress(dataIn),
                    containerType<T>::getSize(dataIn),
                    containerType<T>::getMPItype(), MPI_MAX, MPI_COMM_WORLD);
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
  int numTasks = containerType<V>::getSize(dataOut);
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

#endif
