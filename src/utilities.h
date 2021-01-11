#ifndef UTILS_H
#define UTILS_H

#include <tuple>
#include <string>

// returns the remainder of the division between two integers.
// Note: c++ defines the % operator such that: (a/b)*b + a%b == a    (for b!=0)
// this works as mod() in fortran or % in python when a and b are positive
// But the behavior is not the same for negative integers
// the function below can be used instead
int mod(const int &a, const int &b);

// returns -1 if val<0, 0 if val=0, and 1 if val>0
template<typename T> int sgn(T &val) {
    return (T(0) < val) - (val < T(0));
}

int compress3Indices(const int &i1, const int &i2, const int &i3,
        const int &size1, const int &size2, const int &size3);

std::tuple<int, int, int> decompress3Indices(const int &iTot,
        const int &size1, const int &size2, const int &size3);

int compress2Indices(const int &i1, const int &i2, const int &size1,
        const int &size2);

std::tuple<int, int> decompress2Indices(const int &iTot, const int &size1,
        const int &size2);

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
using CartIndex = NamedType<int, struct CartTag>;

/** BteIndex is used to label the Bloch states entering the BTE
 */
using BteIndex = NamedType<int, struct BteTag>;

/** CalcIndex is used to label the pairs of (chemical potential , temperature)
 * used in transport calculations
 */
using TempIndex = NamedType<int, struct TempTag>;
using ChemPotIndex = NamedType<int, struct ChemPotTag>;

/** WavevectorIndex is used to label k/q points
 */
using WavevectorIndex = NamedType<int, struct WavevectorTag>;

/** BandIndex is used to label bands at a given k/q points
 */
using BandIndex = NamedType<int, struct BandTag>;

/** StateIndex is used to label the Bloch states in the band structure.
 * Not always equal to BteIndex, because some states may be discarded!
 */
using StateIndex = NamedType<int, struct StateTag>;

#endif
