#include "exceptions.h"
#include <string>
#include "periodic_table.h"

int PeriodicTable::findElementByStr(const std::string &speciesName) {
    for (const Element &element : elements) {
        if (element.symbol == speciesName) {
            return element.atomicNumber - 1;
        }
    }
    Error e("Couldn't find element in periodic table");
    return -1;
}

void PeriodicTable::PeriodicTable::setMass(const std::string &speciesName,
        double &x) {
    int i = findElementByStr(speciesName);
    elements[i].mass = x;
}

double PeriodicTable::PeriodicTable::getMass(const std::string &speciesName) {
    int i = findElementByStr(speciesName);
    return elements[i].mass;
}

void PeriodicTable::PeriodicTable::setMassVariance(
        const std::string &speciesName, double &x) {
    int i = findElementByStr(speciesName);
    elements[i].massVariance = x;
}

double PeriodicTable::PeriodicTable::getMassVariance(
        const std::string &speciesName) {
    int i = findElementByStr(speciesName);
    return elements[i].massVariance;
}

int PeriodicTable::getIonicCharge(const std::string &speciesName) {
    int i = findElementByStr(speciesName);
    return elements[i].atomicNumber;
}
