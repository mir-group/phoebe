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
	IO(std::string inputFileName_="", std::string outputFileName_="");
	~IO();
};
