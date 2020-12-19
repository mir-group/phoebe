#include "utilities.h"

long mod(const long &a, const long &b) {
    return (a % b + b) % b;
}

int mod(const int &a, const int &b) {
    return (a % b + b) % b;
}

long compress3Indices(const long &i1, const long &i2, const long &i3,
        const long &size1, const long &size2, const long &size3) {
    (void) size1;
    return i1 * size2 * size3 + i2 * size3 + i3;
}

std::tuple<long, long, long> decompress3Indices(const long &iTot,
        const long &size1, const long &size2, const long &size3) {
    (void) size1;
    long i1 = iTot / (size2 * size3);
    long remainder = iTot - i1 * size2 * size3;
    long i2 = remainder / size3;
    remainder -= i2 * size3;
    long i3 = remainder;
    return {i1,i2,i3};
}

long compress2Indices(const long &i1, const long &i2, const long &size1,
        const long &size2) {
    (void) size1;
    return i1 * size2 + i2;
}

std::tuple<long, long> decompress2Indices(const long &iTot, const long &size1,
        const long &size2) {
    (void) size1;
    long i1 = iTot / size2;
    long i2 = iTot - i1 * size2;
    return {i1,i2};
}


