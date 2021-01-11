#include "utilities.h"

int mod(const int &a, const int &b) {
    return (a % b + b) % b;
}

int compress3Indices(const int &i1, const int &i2, const int &i3,
        const int &size1, const int &size2, const int &size3) {
    (void) size1;
    return i1 * size2 * size3 + i2 * size3 + i3;
}

std::tuple<int, int, int> decompress3Indices(const int &iTot,
        const int &size1, const int &size2, const int &size3) {
    (void) size1;
    int i1 = iTot / (size2 * size3);
    int remainder = iTot - i1 * size2 * size3;
    int i2 = remainder / size3;
    remainder -= i2 * size3;
    int i3 = remainder;
    return {i1,i2,i3};
}

int compress2Indices(const int &i1, const int &i2, const int &size1,
        const int &size2) {
    (void) size1;
    return i1 * size2 + i2;
}

std::tuple<int, int> decompress2Indices(const int &iTot, const int &size1,
        const int &size2) {
    (void) size1;
    int i1 = iTot / size2;
    int i2 = iTot - i1 * size2;
    return {i1,i2};
}


