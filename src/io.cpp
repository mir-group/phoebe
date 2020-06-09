#include "io.h"
#include <algorithm>
#include <exceptions.h>
#include <iomanip>
#include <math.h>
#include <time.h>

char *getCmdOption(char **begin, char **end, const std::string &option) {
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }
  return 0;
}

bool cmdOptionExists(char **begin, char **end, const std::string &option) {
  return std::find(begin, end, option) != end;
}

IO::IO(int argc, char *argv[]) {

  char *inputFileName_ = getCmdOption(argv, argv + argc, "-in");
  char *outputFileName_ = getCmdOption(argv, argv + argc, "-out");
  if (inputFileName_ == nullptr) {
    Error e("Must provide an input file name on command line", 1);
  }

  inputFileName = inputFileName_;
  if (outputFileName_ != nullptr) {
    outputFileName = outputFileName_;
    outputFile.open(outputFileName);
    std::streambuf *coutbuf = std::cout.rdbuf(); // save old buf
    (void)coutbuf;                       // suppress unused variable error
    std::cout.rdbuf(outputFile.rdbuf()); // redirect std::cout to out.txt!
  }
};

IO::~IO() {
  if (outputFile.is_open()) {
    outputFile.close();
  }
}

std::string IO::getInputFileName() { return inputFileName; }

void IO::welcome() {
  std::string welcomeMsg =
      "\n"
      "  8888888b.  888                        888\n"
      "  888   Y88b 888                        888\n"
      "  888    888 888                        888\n"
      "  888   d88P 88888b.   .d88b.   .d88b.  88888b.   .d88b.  \n"
      "  8888888P'  888 '88b d88''88b d8P  Y8b 888 '88b d8P  Y8b \n"
      "  888        888  888 888  888 88888888 888  888 88888888 \n"
      "  888        888  888 Y88..88P Y8b.     888 d88P Y8b.     \n"
      "  888        888  888  'Y88P'   'Y8888  88888P'   'Y8888' \n"
      "\n";
  std::cout << welcomeMsg;
}

void IO::goodbye() { std::cout << "Exiting program" << std::endl; }

LoopPrint::LoopPrint(const std::string &task_, const std::string step_,
                     const long &numSteps_) {
  task = task_;
  step = step_;
  numSteps = numSteps_;

  long numRep = 10; // number of intermediate reports
  if (numSteps < numRep) {
    reportEvery = 1;
  } else {
    reportEvery = numSteps / numRep;
  }

  // note: time returns the time (in secs) elapsed from 1st Jan 1970
  initialTime = std::chrono::steady_clock::now();

  std::cout << "\n";
  std::cout << "Started " << task << " with " << numSteps << " " << step
            << ".\n";

  stepDigits = long(log10(numSteps)) + 1; // number of digits in numSteps
}

void LoopPrint::update() {
  // Update prediction info for current task.

  currentStep += 1;

  if (currentStep <= 2 || (currentStep + 1) % reportEvery == 0 ||
      currentStep == numSteps - 1) {

    time_point currentTime;
    currentTime = std::chrono::steady_clock::now();

    // format currentTime nicely
    char s[200];
    time_t curr = time(NULL);
    struct tm *p = localtime(&curr);
    strftime(s, 200, "%F, %T", p);

    time_delta elapsedTime = currentTime - initialTime;
    double timeLeft;

    if (currentStep == 2) { // we compare with the third step
      deltaTime = elapsedTime;
    } else if (currentStep > 2) {
      timeLeft = std::chrono::duration_cast<std::chrono::nanoseconds>(
                     elapsedTime - deltaTime)
                     .count() /
                 1e9 / (currentStep - 2.) * (numSteps - currentStep + 1.);
    }

    if ((currentStep == 0 || currentStep == 2 || currentStep == numSteps - 1) ||
        (currentStep + 1) % reportEvery == 0) {

      long percentage = double(currentStep + 1) / numSteps * 100.;

      std::cout << s << " | ";
      std::cout << std::setw(3) << percentage << "% | ";
      std::cout << std::setw(stepDigits) << currentStep + 1
                << std::setw(stepDigits) << " / " << numSteps;
      if (currentStep > 2) {
        std::cout << " | remaining: " << std::setw(8) << std::setprecision(2)
                  << std::scientific << timeLeft << std::fixed << " s.\n";
      } else {
        std::cout << "\n";
      }
    }
    std::cout << std::resetiosflags(std::cout.flags());
  }
}

void LoopPrint::close() {
  // print timing results
  time_point currentTime;
  currentTime = std::chrono::steady_clock::now();
  std::cout << "Elapsed time: " << std::setprecision(3)
            << std::chrono::duration_cast<std::chrono::nanoseconds>(
                   currentTime - initialTime)
                       .count() /
                   1e9
            << " s.\n";
}
