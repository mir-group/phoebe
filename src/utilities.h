#ifndef UTILS_H
#define UTILS_H

#include <tuple>
#include <string>

// returns the remainder of the division between two integers.
// Note: c++ defines the % operator such that: (a/b)*b + a%b == a    (for b!=0)
// this works as mod() in fortran or % in python when a and b are positive
// But the behavior is not the same for negative integers
// the function below can be used instead
long mod(const long &a, const long &b);
int mod(const int &a, const int &b);

// checks if string ends with a suffix
bool hasSuffix(const std::string &str, const std::string &suffix);

// returns -1 if val<0, 0 if val=0, and 1 if val>0
template<typename T> int sgn(T &val) {
    return (T(0) < val) - (val < T(0));
}

long compress3Indices(const long &i1, const long &i2, const long &i3,
        const long &size1, const long &size2, const long &size3);

std::tuple<long, long, long> decompress3Indices(const long &iTot,
        const long &size1, const long &size2, const long &size3);

long compress2Indices(const long &i1, const long &i2, const long &size1,
        const long &size2);

std::tuple<long, long> decompress2Indices(const long &iTot, const long &size1,
        const long &size2);

// A function to allocate a dynamically sized array. It tricks the 
// compiler into thinking the size is a constant via the const identifier
// on the argument. This resolves issues with VLAs -- see crystal.cpp
template <typename T> T* allocate(T *&array, const unsigned int size){
        array = new T [size];
        return array;
}

/** Class for implementing strong typing.
 * In fact, when using the methods (de)compress2(3)indices, it's easy to
 * mix the order of indices. This class can be used to enforce the type
 * of the index, without incurring in a penalty in the code speed.
 * Taken from https://www.fluentcpp.com/2017/05/05/news-strong-types-are-free/
 */
template<typename T, typename Parameter>
class NamedType {
public:
    explicit NamedType(T const &value) :
            value_(value) {
    }
    T& get() {
        return value_;
    }
    T const& get() const {
        return value_;
    }
private:
    T value_;
};

/** CartIndex is used to label cartesian directions
 */
using CartIndex = NamedType<long, struct CartTag>;

/** BteIndex is used to label the Bloch states entering the BTE
 */
using BteIndex = NamedType<long, struct BteTag>;

/** CalcIndex is used to label the pairs of (chemical potential , temperature)
 * used in transport calculations
 */
using CalcIndex = NamedType<long, struct CalcTag>;
using TempIndex = NamedType<long, struct TempTag>;
using ChemPotIndex = NamedType<long, struct ChemPotTag>;

/** WavevectorIndex is used to label k/q points
 */
using WavevectorIndex = NamedType<long, struct WavevectorTag>;

/** BandIndex is used to label bands at a given k/q points
 */
using BandIndex = NamedType<long, struct BandTag>;

/** StateIndex is used to label the Bloch states in the bandstructure.
 * Not always equal to BteIndex, because some states may be discarded!
 */
using StateIndex = NamedType<long, struct StateTag>;

#endif