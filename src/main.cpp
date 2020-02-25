#include <iostream>
#include "periodic_table.h"

//using namespace periodicTable;

int main() {
	using namespace pt;

	std::cout << pt::periodicTable[7].symbol + "\n";
	std::cout << pt::periodicTable[7].mass;
	std::cout << "\n";
	std::cout << pt::periodicTable[7].massVariance;
	std::cout << "\n";

	return(0);
}
