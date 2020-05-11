#include <iterator>
#include <sstream>
#include <vector>

long mod(long a, long b) {
	return ( a%b + b ) % b;
}

bool hasSuffix(const std::string & str, const std::string & suffix) {
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

long compressIndeces(long i, long j, long k, long size1, long size2,
		long size3) {
	(void) size1;
	return i*size2*size3 + j*size3 + k;
}
