#include "io.h"
#include <algorithm>
#include <exceptions.h>
#include <time.h>
#include <iomanip>
#include <math.h>
#include "io.h"
#include "mpiHelper.h"

// Utility to get the command line option from it's name
char* getCmdOption(char **begin, char **end, const std::string &option) {
    char **itr = std::find(begin, end, option);
    if (itr != end && ++itr != end) {
        return *itr;
    }
    return 0;
}

// utility to check if the command line option exists
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
    //redirect std::cout to outputFilename, if passed on command line
    if (outputFileName_ != nullptr) {
        outputFileName = outputFileName_;
        outputFile.open(outputFileName);
        std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
        (void) coutbuf; // suppress unused variable error
        std::cout.rdbuf(outputFile.rdbuf()); //redirect std::cout to outputFile
    }
}
;

IO::~IO() {
    if (outputFile.is_open()) {
        outputFile.close();
    }
}

std::string IO::getInputFileName() {
    return inputFileName;
}

void IO::welcome() {
    if ( ! mpi->mpiHead() ) return;
    std::string welcomeMsg = "\n"
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

void IO::goodbye() {
    if ( ! mpi->mpiHead() ) return;
    std::cout << "Exiting program" << std::endl;
}

LoopPrint::LoopPrint(const std::string &task_, const std::string step_,
        const long &numSteps_) {
    if ( ! mpi->mpiHead() ) return;

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
    initialTime = time(NULL);

    std::cout << "\n";
    std::cout << "Started " << task << " with " << numSteps << " " << step
            << ".\n";

    stepDigits = long(log10(numSteps)) + 1; // number of digits in numSteps
}

void LoopPrint::update() {
    // Update prediction info for current task.
    if ( ! mpi->mpiHead() ) return;

    currentStep += 1;

    // we report at the first three steps, the last step, and ~10% intervals
    if (currentStep <= 2 || (currentStep + 1) % reportEvery == 0
            || currentStep == numSteps - 1) {

        // get system time
        time_t currentTime;
        currentTime = time(NULL);

        // format currentTime nicely
        char s[200];
        struct tm *p = localtime(&currentTime);
        strftime(s, 200, "%F, %T", p);

        // time from the beginning of LoopPrint
        time_t elapsedTime = currentTime - initialTime;
        time_t timeLeft;

        if (currentStep == 2) { // we compare with the third step
            // the first two steps are used to make an estimate later
            deltaTime = elapsedTime;
        } else if (currentStep > 2) { // make the estimate
            // estimate the remaining time
            timeLeft = (elapsedTime - deltaTime) / (currentStep - 2.)
                    * (numSteps - currentStep + 1.);
        }

        // print to screen the loop info
        if ((currentStep == 0 || currentStep == 2 || currentStep == numSteps - 1)
                || (currentStep + 1) % reportEvery == 0) {

            long percentage = double(currentStep + 1) / numSteps * 100.;

            std::cout << s << " | ";
            std::cout << std::setw(3) << percentage << "% | ";
            std::cout << std::setw(stepDigits) << currentStep + 1
                    << std::setw(stepDigits) << " / " << numSteps;
            if (currentStep > 2) {
                std::cout << " | remaining: " << timeLeft << " s.\n";
            } else {
                std::cout << "\n";
            }
        }
    }
}

void LoopPrint::close() {
    if ( ! mpi->mpiHead() ) return;
    // print timing results
    time_t currentTime;
    currentTime = time(NULL);
    std::cout << "Elapsed time: " << currentTime - initialTime << " s.\n";
}

