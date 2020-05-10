#include <string>
#include <fstream>
#include <cstdarg>
#include <iostream>
#include <cstdio>

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
