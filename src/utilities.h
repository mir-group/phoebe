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

/** Class for implementing strong typing.
 * In fact, when using the methods (de)compress2(3)indices, it's easy to
 * mix the order of indices. This class can be used to enforce the type
 * of the index, without incurring in a penalty in the code speed.
 * Taken from https://www.fluentcpp.com/2017/05/05/news-strong-types-are-free/
 */
template <typename T, typename Parameter>
class NamedType
{
public:
    explicit NamedType(T const& value) : value_(value) {}
    T& get() { return value_; }
    T const& get() const {return value_; }
private:
    T value_;
};

using DimIndex = NamedType<double, struct DimTag>;
using CalcIndex = NamedType<double, struct CalcTag>;
using TempIndex = NamedType<double, struct TempTag>;
using ChemPotIndex = NamedType<double, struct ChemPotTag>;
using WavevectorIndex = NamedType<double, struct WavevectorTag>;
using BandIndex = NamedType<double, struct BandTag>;

#endif
