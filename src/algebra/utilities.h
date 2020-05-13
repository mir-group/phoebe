#ifndef UTILS_H
#define UTILS_H

#include <tuple>
#include <string>

// returns the remainder of the division between two integers.
// Note: c++ defines the % operator such that: (a/b)*b + a%b == a    (for b!=0)
// this works as mod() in fortran or % in python when a and b are positive
// But the behavior is not the same for negative integers
// the function below can be used instead
long mod(const long & a, const long & b);

// checks if string ends with a suffix
bool hasSuffix(const std::string & str, const std::string & suffix);

// returns -1 if val<0, 0 if val=0, and 1 if val>0
template <typename T> int sgn(T & val) {
    return (T(0) < val) - (val < T(0));
}

long compress3Indeces(const long & i1, const long & i2, const long & i3,
		const long & size1, const long & size2, const long & size3);

std::tuple<long,long,long> decompress3Indeces(const long & iTot,
		const long & size1, const long & size2, const long & size3);

long compress2Indeces(const long & i1, const long & i2, const long & size1,
		const long & size2);

std::tuple<long,long> decompress2Indeces(const long & iTot, const long & size1,
		const long & size2);


#endif
