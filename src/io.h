#ifndef IO_H
#define IO_H

#include <chrono>
#include <cstdarg>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>
#include <time.h>
#include <type_traits>

class IO {
private:
  std::ofstream outputFile;
  std::string outputFileName = "";
  std::string inputFileName = "";

public:
  IO(int argc, char *argv[]);
  ~IO();
  std::string getInputFileName();
  void welcome();
  void goodbye();
};

// object to report time progress in a loop
class LoopPrint {
private:
public:
  LoopPrint(const std::string &task, const std::string step,
            const long &numSteps);
  void update();
  void close();

private:
  typedef std::chrono::steady_clock::time_point time_point;
  typedef std::chrono::steady_clock::duration time_delta;
  long reportEvery;
  long numSteps;
  std::string task;
  std::string step;
  long currentStep = -1;
  time_point initialTime;
  time_delta deltaTime;
  long stepDigits;
};

#endif
