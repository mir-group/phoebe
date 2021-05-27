#include "utilities.h"
#include <fstream>
#include <iostream>
#include <mpiHelper.h>
#include <unistd.h>
#include <iomanip>

int mod(const int &a, const int &b) { return (a % b + b) % b; }

int compress3Indices(const int &i1, const int &i2, const int &i3,
                     const int &size1, const int &size2, const int &size3) {
  (void)size1;
  return i1 * size2 * size3 + i2 * size3 + i3;
}

std::tuple<int, int, int> decompress3Indices(const int &iTot, const int &size1,
                                             const int &size2,
                                             const int &size3) {
  (void)size1;
  int i1 = iTot / (size2 * size3);
  int remainder = iTot - i1 * size2 * size3;
  int i2 = remainder / size3;
  remainder -= i2 * size3;
  int i3 = remainder;
  return {i1, i2, i3};
}

int compress2Indices(const int &i1, const int &i2, const int &size1,
                     const int &size2) {
  (void)size1;
  return i1 * size2 + i2;
}

std::tuple<int, int> decompress2Indices(const int &iTot, const int &size1,
                                        const int &size2) {
  (void)size1;
  int i1 = iTot / size2;
  int i2 = iTot - i1 * size2;
  return {i1, i2};
}

std::tuple<double, double> memoryUsage() {
  double vm_usage = 0.;
  double resident_set = 0.;
  std::ifstream stat_stream("/proc/self/stat",
                            std::ios_base::in); // get info from proc directory
  // create some variables to get info
  std::string pid, comm, state, ppid, pgrp, session, tty_nr;
  std::string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  std::string utime, stime, cutime, cstime, priority, nice;
  std::string O, itrealvalue, starttime;
  unsigned long vsize;
  long rss;
  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
      tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
      stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
      starttime >> vsize >> rss; // don't care about the rest
  stat_stream.close();
  long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                      1024; // for x86-64 is configured to use 2MB pages
  vm_usage = vsize / 1024.0;
  resident_set = rss * page_size_kb;

  mpi->allReduceSum(&vm_usage);
  mpi->allReduceSum(&resident_set);
  if ( mpi->mpiHead() ) {
    std::cout << std::setprecision(4);
    std::cout << "Snapshot of Phoebe's memory usage:\n";
    std::cout << "VM: " << vm_usage/1024./1024. << " (GB). RSS: "
              << resident_set/1024./1024. << " (GB)\n" << std::endl;
  }

  return {vm_usage, resident_set};
}
