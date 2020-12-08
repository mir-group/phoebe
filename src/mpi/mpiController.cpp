#include "mpiController.h"
#include <chrono>
#include <iostream>
#include <vector>

#ifdef MPI_AVAIL
#include <mpi.h>
#endif

// default constructor
MPIcontroller::MPIcontroller() {

#ifdef MPI_AVAIL
  // start the MPI environment
  int errCode = MPI_Init(NULL, NULL);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }

  // set this so that MPI returns errors and lets us handle them, rather
  // than using the default, MPI_ERRORS_ARE_FATAL
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN);

  // get the number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // get rank of current process
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // start a timer
  startTime = MPI_Wtime();

#else
  // To maintain consistency when running in serial
  size = 1;
  rank = 0;
  startTime = std::chrono::steady_clock::now();
#endif
}

void MPIcontroller::finalize() const {
#ifdef MPI_AVAIL
  barrier();
  if (mpiHead()) {
    fprintf(stdout, "Final time: %3f\n ", MPI_Wtime() - startTime);
  }
  MPI_Finalize();
#else
  std::cout << "Final time: "
            << std::chrono::duration_cast<std::chrono::microseconds>(
                   std::chrono::steady_clock::now() - startTime)
                       .count() *
                   1e-6
            << " secs" << std::endl;
#endif
}

// Utility functions  -----------------------------------

// get the error string and print it to stderr before returning
void MPIcontroller::errorReport(int errCode) const {
#ifdef MPI_AVAIL
  char errString[BUFSIZ];
  int lengthOfString;

  MPI_Error_string(errCode, errString, &lengthOfString);
  fprintf(stderr, "Error from rank %3d: %s\n", rank, errString);
  MPI_Abort(MPI_COMM_WORLD, errCode);
#endif
}

void MPIcontroller::time() const {
#ifdef MPI_AVAIL
  fprintf(stdout, "Time for rank %3d : %3f\n", rank, MPI_Wtime() - startTime);
#else
  std::cout << "Time for rank 0 :"
            << std::chrono::duration_cast<std::chrono::seconds>(
                   std::chrono::steady_clock::now() - startTime)
                   .count()
            << " secs" << std::endl;
#endif
}
// Asynchronous support functions -----------------------------------------
void MPIcontroller::barrier() const {
#ifdef MPI_AVAIL
  int errCode;
  errCode = MPI_Barrier(MPI_COMM_WORLD);
  if (errCode != MPI_SUCCESS) {
    errorReport(errCode);
  }
#endif
}

// Labor division functions -----------------------------------------
std::vector<long> MPIcontroller::divideWork(size_t numTasks) {
  // return a vector of the start and stop points for task division
  std::vector<long> divs(2);
  divs[0] = (numTasks * rank) / size;
  divs[1] = (numTasks * (rank + 1)) / size;
  return divs;
}

std::vector<long> MPIcontroller::divideWorkIter(size_t numTasks) {
  // return a vector of indices for tasks to be completed by thsi process 
  int start = (numTasks * rank) / size;
  int stop = (numTasks * (rank + 1)) / size;
  size_t localSize = stop - start;
  std::vector<long> divs(localSize);
  for (long i = start; i < stop; i++) {
    divs[i - start] = i;
  }
  return divs;
}

// Helper function to re-establish work divisions for MPI calls requiring
// the number of tasks given to each point
std::tuple<std::vector<int>, std::vector<int>>
MPIcontroller::workDivHelper(size_t numTasks) const {

  std::vector<int> workDivs(size);
  // start points for each rank's work
  std::vector<int> workDivisionHeads(size);
  std::vector<int> workDivisionTails(size);
  // Recreate work division instructions
  for (int i = 0; i < size; i++) {
    workDivisionHeads[i] = (numTasks * i) / size;
    workDivisionTails[i] = (numTasks * (i + 1)) / size;
  }
  /** Note: it is important to compute workDivs as the subtraction of two
   * other variables. Some compilers (e.g. gcc 9.3.0 on Ubuntu) may optimize
   * the calculation of workDivs setting it to workDivs[i]=numTasks/size ,
   * which doesn't work when the division has a remainder.
   */
  for (int i = 0; i < size; i++) {
    workDivs[i] = workDivisionTails[i] - workDivisionHeads[i];
  }
  return {workDivs, workDivisionHeads};
}
