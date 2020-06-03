#include <string>
#include <fstream>
#include <cstdarg>
#include <iostream>
#include <cstdio>
#include <time.h>

class IO {
private:
	std::ofstream outputFile;
	std::string outputFileName = "";
	std::string inputFileName = "";
public:
	IO(int argc, char* argv[]);
	~IO();
	std::string getInputFileName();
	void welcome();
	void goodbye();
};

// object to report time progress in a loop
class LoopPrint {
private:
public:
	LoopPrint(const std::string & task, const std::string step,
			const long & numSteps);
	void update();
	void close();
private:
	long reportEvery;
	long numSteps;
	std::string task;
	std::string step;
	long currentStep = -1;
	time_t initialTime;
	long deltaTime;
	long stepDigits;
};
