#include "context.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include "constants.h"
#include "exceptions.h"

/** Returns true if a string contains a substring.
 */
bool patternInString(const std::string &s, const std::string &pattern) {
  if (s.find(pattern) != std::string::npos) {
    return true;
  } else {
    return false;
  }
}

/** Parse a string of format "key = value" to return a boolean value.
 */
bool parseBool(std::string line) {
  std::string delimeter = "=";
  size_t pos = line.find(delimeter);
  std::string s = line.substr(pos + 1);

  std::transform(s.begin(), s.end(), s.begin(), ::tolower);
  s.erase(
      std::remove_if(s.begin(), s.end(), [](char c) { return !isalpha(c); }),
      s.end());

  if (s == "false") {
    return false;
  } else if (s == "true") {
    return true;
  } else if (s == "0") {
    return false;
  } else if (s == "1") {
    return true;
  } else {
    Error e("Couldn't fix boolean value while parsing");
    return false;
  }
};

/** Parse a string of format "key = value" to return a double value.
 */
double parseDouble(std::string line) {
  std::string delimeter = "=";
  size_t pos = line.find(delimeter);
  std::string value = line.substr(pos + 1);
  return std::stod(value); // convert to double
};

/** Parse a string of format "key = value units" to return a double value
 * converted in rydberg atomic units.
 */
double parseDoubleWithUnits(std::string line) {
  double x;

  std::string delimeter = "=";
  size_t pos = line.find(delimeter);
  std::string value = line.substr(pos + 1);
  x = std::stod(value); // convert to double

  // now check the units and convert
  if (patternInString(line, "eV")) {
    x /= energyRyToEv;
  }
  if (patternInString(line, "cmm1")) {
    x /= ryToCmm1;
  }
  if (patternInString(line, "ps")) {
    x /= timeRyToFs * 1.0e-3;
  }
  if (patternInString(line, "fs")) {
    x /= timeRyToFs;
  }
  if (patternInString(line, "mum")) {
    x /= distanceBohrToMum;
  }

  return x;
};

/** Parse a string of format "key = [value1,value2]" to return a vector double.
 */
std::vector<double> parseDoubleList(std::string line) {
  std::string delimeter = "[";
  size_t pos1 = line.find_first_of(delimeter);
  delimeter = "]";
  size_t pos2 = line.find_last_of(delimeter);

  if (pos1 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }
  if (pos2 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }

  std::string s = line.substr(pos1 + 1, pos2 - pos1 - 1);

  std::vector<double> x;
  delimeter = ",";
  while ((pos1 = s.find(delimeter)) != std::string::npos) {
    std::string token = s.substr(0, pos1);
    x.push_back(std::stod(token));
    s.erase(0, pos1 + delimeter.length());
  }
  // Must not forget the last element in the list
  x.push_back(std::stod(s));

  return x;
};

/** Parse a string of format "key = [value1,value2]" to return a vector double.
 */
std::vector<int> parseIntList(std::string line) {
  std::string delimeter = "[";
  size_t pos1 = line.find_first_of(delimeter);
  delimeter = "]";
  size_t pos2 = line.find_last_of(delimeter);

  if (pos1 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }
  if (pos2 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }

  std::string s = line.substr(pos1 + 1, pos2 - pos1 - 1);

  std::vector<int> x;
  delimeter = ",";
  while ((pos1 = s.find(delimeter)) != std::string::npos) {
    std::string token = s.substr(0, pos1);
    x.push_back(std::stoi(token));
    s.erase(0, pos1 + delimeter.length());
  }
  // Must not forget the last element in the list
  x.push_back(std::stoi(s));

  return x;
};

/** Parse a string of format "key = value units" to return an integer value.
 */
long parseLong(std::string line) {
  std::string delimeter = "=";
  size_t pos = line.find(delimeter);
  std::string value = line.substr(pos + 1);
  return std::stoi(value); // convert to integer
};

/** Parse a string of format "key = [val1,val2]" to return a vector of ints.
 */
std::vector<long> parseLongList(std::string line) {
  std::string delimeter = "[";
  size_t pos1 = line.find_first_of(delimeter);
  delimeter = "]";
  size_t pos2 = line.find_last_of(delimeter);

  if (pos1 == std::string::npos) {
    Error e("Error in parseLongList");
  }
  if (pos2 == std::string::npos) {
    Error e("Error in parseLongList");
  }

  std::string s = line.substr(pos1 + 1, pos2 - pos1 - 1);
  delimeter = ",";
  std::vector<long> x;
  while ((pos1 = s.find(delimeter)) != std::string::npos) {
    std::string token = s.substr(0, pos1);
    x.push_back(std::stoi(token)); // convert to integer
    s.erase(0, pos1 + delimeter.length());
  }
  // Must not forget the last element in the list
  x.push_back(std::stoi(s));

  return x;
};

/** Parse a string of format "key = value" to return a string value.
 */
std::string parseString(std::string line) {
  std::string delimeter = "'";
  size_t pos1 = line.find_first_of(delimeter);
  size_t pos2 = line.find_last_of(delimeter);

  if (pos1 == std::string::npos) {
    delimeter = "\"";
    pos1 = line.find_first_of(delimeter);
    pos2 = line.find_last_of(delimeter);
    if (pos1 == std::string::npos) {
      Error e("Couldn't solve string parsing");
    }
  }

  if (pos1 == pos2) {
    Error e("Error parsing string from user input");
  }
  std::string x = line.substr(pos1 + 1, pos2 - pos1 - 1);
  return x;
};

/** Parse a string of format "key = [val1,val2]" to return a vector of strings.
 */
std::vector<std::string> parseStringList(std::string line) {
  // remove empty spaces
  line.erase(std::remove(line.begin(), line.end(), ' '), line.end());

  std::string delimeter;
  delimeter = "[";
  size_t pos1 = line.find_first_of(delimeter);
  delimeter = "]";
  size_t pos2 = line.find_last_of(delimeter);

  if (pos1 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }
  if (pos2 == std::string::npos) {
    Error e("Error in parseDoubleList");
  }

  std::string s = line.substr(pos1 + 1, pos2 - pos1 - 1);
  delimeter = ",";
  std::vector<std::string> x;
  while ((pos1 = s.find(delimeter)) != std::string::npos) {
    std::string token = s.substr(0, pos1);
    // we also remove the " symbols
    token.erase(std::remove(token.begin(), token.end(), '"'), token.end());
    x.push_back(token);
    s.erase(0, pos1 + delimeter.length());
  }
  // Must not forget the last element in the list
  s.erase(std::remove(s.begin(), s.end(), '"'), s.end());
  x.push_back(s);

  return x;
};

/** Parse the block of information on the crystal structure.
 * Format:
 * Atom1Name   cartesianCoordx   cartesianCoordy   cartesianCoordz
 * Atom2Name   cartesianCoordx   cartesianCoordy   cartesianCoordz
 * ...
 * The coordinates must be provided in Angstroms.
 */
std::tuple<Eigen::MatrixXd, Eigen::VectorXi, std::vector<std::string>>
parseCrystal(std::vector<std::string> &lines) {
  long numAtoms = lines.size();
  Eigen::MatrixXd atomicPositions(numAtoms, 3);
  Eigen::VectorXi atomicSpecies(numAtoms);
  std::vector<std::string> speciesNames;

  int counter = 0;
  for (std::string line : lines) {
    // split line by spaces
    std::stringstream ss(line);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> splitLine(begin, end);

    std::string thisElement = splitLine[0];
    // if new element, insert into list of species
    if (std::find(speciesNames.begin(), speciesNames.end(), thisElement) ==
        speciesNames.end()) {
      // thisElement not in speciesNames, add it
      speciesNames.push_back(thisElement);
    }
    // find the index of the current element
    long index = 0;
    for (auto speciesName : speciesNames) {
      if (speciesName == thisElement) {
        break;
      }
      index += 1;
    }
    // save species and positions
    atomicSpecies(counter) = index;
    atomicPositions(counter, 0) = std::stod(splitLine[1]);
    atomicPositions(counter, 1) = std::stod(splitLine[2]);
    atomicPositions(counter, 2) = std::stod(splitLine[3]);

    counter++;
  }
  atomicPositions /= distanceBohrToAng;
  return {atomicPositions, atomicSpecies, speciesNames};
}

/** Reads the block with the information on the path of points in the
 * Brillouin zone.
 *
 * Format:
 * """
 * PointName1  coord1X  coord1Y  coord1Z  PointName2  coord2X  coord2Y  coord2Z
 * PointName2  coord2X  coord2Y  coord2Z  PointName3  coord3X  coord3Y  coord3Z
 * """
 *
 * Each line identifies a segment of the path in the Brillouin zone.
 * PointName1 and PointName2 are the names of the two special points at the
 * opposite ends of the segment.
 * coord1/2 are the crystal coordinates of these special points.
 * The code will pad the segments with equidistant points, as specified by the
 * parameter deltaPath (the distance between points).
 */
std::tuple<std::vector<std::string>, Eigen::Tensor<double, 3>>
parsePathExtrema(std::vector<std::string> &lines) {

  long numSegments = lines.size();
  Eigen::Tensor<double, 3> pathExtrema(numSegments, 2, 3);
  pathExtrema.setZero();
  std::vector<std::string> pathLabels;

  long i = 0;
  for (std::string line : lines) {
    // split line by spaces
    std::stringstream ss(line);
    std::istream_iterator<std::string> begin(ss);
    std::istream_iterator<std::string> end;
    std::vector<std::string> splitLine(begin, end);

    pathLabels.push_back(splitLine[0]);
    pathLabels.push_back(splitLine[4]);

    pathExtrema(i, 0, 0) = std::stod(splitLine[1]);
    pathExtrema(i, 0, 1) = std::stod(splitLine[2]);
    pathExtrema(i, 0, 2) = std::stod(splitLine[3]);

    pathExtrema(i, 1, 0) = std::stod(splitLine[5]);
    pathExtrema(i, 1, 1) = std::stod(splitLine[6]);
    pathExtrema(i, 1, 2) = std::stod(splitLine[7]);

    i++;
  }
  return {pathLabels, pathExtrema};
}

/** Parse an input block (e.g. for crystal or points path).
 * Returns empty strings/vectors if no block name is found at line "lineCounter"
 * of the input file.
 * Otherwise, returns the block name (i.e. what's named after "begin")
 * and the lines inside the block.
 */
std::tuple<std::string, std::vector<std::string>>
parseBlockNameValue(const std::vector<std::string> &lines,
                    const int &lineCounter) {
  std::string line;
  line = lines[lineCounter];
  if (!patternInString(line, "begin")) {
    std::string empty1 = "";
    std::vector<std::string> empty2;
    return {empty1, empty2};
  } else {
    std::string pattern = "begin";
    std::size_t found = line.find(pattern);
    std::string blockName = line.substr(found + pattern.size() + 1);
    std::vector<std::string> val;
    for (int unsigned i = lineCounter + 1; i < lines.size(); i++) {
      if (patternInString(lines[i], "end")) {
        break;
      }
      val.push_back(lines[i]);
    }
    return {blockName, val};
  }
}

std::vector<std::string> &Context::split(const std::string &s, char delim,
                                         std::vector<std::string> &elems) {
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    if (item.length() > 0) {
      elems.push_back(item);
    }
  }
  return elems;
}

/** Split a string by a char delimeter.
 */
std::vector<std::string> Context::split(const std::string &s, char delim) {
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

/** Checks if the line of the input file contains a key=value statement.
 */
bool lineHasParameter(const std::string &line) {
  if (std::find(line.begin(), line.end(), '=') != line.end()) {
    return true;
  } else {
    return false;
  }
}

/** Given a line of the form "key = value", returns the string "key" and the
 * string containing "value". "value" could contain a string, a list or other
 * things to be further parsed.
 */
std::tuple<std::string, std::string>
parseParameterNameValue(const std::string &line) {
  // we assume that there is "=" in the string
  std::string sep = "=";
  size_t position = line.find(sep); // unsigned integer
  std::string s = line.substr(0, position);
  s.erase(std::remove_if(s.begin(), s.end(), ::isspace), s.end());

  std::string val = line.substr(position + 1, line.size());

  return {s, val};
}

void Context::setupFromInput(std::string fileName) {
  std::vector<std::string> lines;
  std::string line;

  // open input file and read content
  std::ifstream infile(fileName);
  while (std::getline(infile, line)) {
    std::vector<std::string> tokens = split(line, ';');
    for (std::string t : tokens) {
      lines.push_back(t);
    }
  }
  infile.close();

  // there are

  int lineCounter = 0;
  for (std::string line : lines) {
    if (line.empty()) { // nothing to do
      continue;

      // line with pair (key,value)
    } else if (lineHasParameter(line)) {
      auto tup = parseParameterNameValue(line);
      auto parameterName = std::get<0>(tup);
      auto val = std::get<1>(tup);

      if (parameterName == "phD2FileName") {
        phD2FileName = parseString(val);
      }

      if (parameterName == "phD3FileName") {
        phD3FileName = parseString(val);
      }

      if (parameterName == "sumRuleD2") {
        sumRuleD2 = parseString(val);
      }

      if (parameterName == "electronH0Name") {
        electronH0Name = parseString(val);
      }

      if (parameterName == "wannier90Prefix") {
        wannier90Prefix = parseString(val);
      }

      if (parameterName == "quantumEspressoPrefix") {
        quantumEspressoPrefix = parseString(val);
      }

      if (parameterName == "elPhInterpolation") {
        elPhInterpolation = parseString(val);
      }

      if (parameterName == "epwFileName") {
        setEpwFileName(parseString(val));
      }

      if (parameterName == "electronFourierCutoff") {
        double x = parseDouble(val);
        electronFourierCutoff = x;
      }

      if (parameterName == "qMesh") {
        std::vector<long> vecMesh = parseLongList(val);
        qMesh(0) = vecMesh[0];
        qMesh(1) = vecMesh[1];
        qMesh(2) = vecMesh[2];
      }

      if (parameterName == "kMesh") {
        std::vector<long> vecMesh = parseLongList(val);
        kMesh(0) = vecMesh[0];
        kMesh(1) = vecMesh[1];
        kMesh(2) = vecMesh[2];
      }

      if (parameterName == "windowType") {
        windowType = parseString(val);
      }

      if (parameterName == "windowEnergyLimit") {
        std::vector<double> winLim = parseDoubleList(val);
        windowEnergyLimit[0] = std::min(winLim[0], winLim[1]) / energyRyToEv;
        windowEnergyLimit[1] = std::max(winLim[0], winLim[1]) / energyRyToEv;
      }

      if (parameterName == "windowPopulationLimit") {
        windowPopulationLimit = parseDouble(val);
      }

      if (parameterName == "chemicalPotentials") {
        std::vector<double> x = parseDoubleList(val);
        chemicalPotentials = Eigen::VectorXd::Zero(x.size());
        for (long unsigned i = 0; i < x.size(); i++) {
          chemicalPotentials(i) = x[i] / energyRyToEv;
        }
      }

      if (parameterName == "dopings") {
        std::vector<double> x = parseDoubleList(val);
        dopings = Eigen::VectorXd::Zero(x.size());
        for (long unsigned i = 0; i < x.size(); i++) {
          dopings(i) = x[i];
        }
      }

      if (parameterName == "temperatures") {
        std::vector<double> x = parseDoubleList(val);
        temperatures = Eigen::VectorXd::Zero(x.size());
        for (long unsigned i = 0; i < x.size(); i++) {
          temperatures(i) = x[i] / temperatureAuToSi;
        }
      }

      if (parameterName == "appName") {
        appName = parseString(val);
      }

      if (parameterName == "solverBTE") {
        solverBTE = parseStringList(val);
      }

      if (parameterName == "convergenceThresholdBTE") {
        convergenceThresholdBTE = parseDouble(val);
      }

      if (parameterName == "maxIterationsBTE") {
        maxIterationsBTE = parseLong(val);
      }

      if (parameterName == "dimensionality") {
        dimensionality = parseLong(val);
      }

      if (parameterName == "dosMinEnergy") {
        dosMinEnergy = parseDoubleWithUnits(val);
      }

      if (parameterName == "dosMaxEnergy") {
        dosMaxEnergy = parseDoubleWithUnits(val);
      }

      if (parameterName == "dosDeltaEnergy") {
        dosDeltaEnergy = parseDoubleWithUnits(val);
      }

      if (parameterName == "deltaPath") {
        deltaPath = parseDouble(val);
      }

      if (parameterName == "fermiLevel") {
        fermiLevel = parseDoubleWithUnits(val);
      }

      if (parameterName == "hasSpinOrbit") {
        hasSpinOrbit = parseBool(val);
      }

      if (parameterName == "numOccupiedStates") {
        // note: numOccupiedStates refers to the number of states that are
        // occupied
        // for Wannier: the number of Wannier states that are full
        // for Fourier: the number of occupied bands
        // remember to NOT count the spin degeneracy
        double x = parseDouble(val);
        if (!hasSpinOrbit)
          x *= 2;
        numOccupiedStates = x;
      }

      if (parameterName == "smearingMethod") {
        std::string x_ = parseString(val);
        // TODO: this is hardcoded, should be fixed how we validate input
        if (x_ == "gaussian") {
          smearingMethod = 0;
        } else if (x_ == "adaptiveGaussian") {
          smearingMethod = 1;
        } else if (x_ == "tetrahedron") {
          smearingMethod = 2;
        } else {
          smearingMethod = -1;
        }
      }

      if (parameterName == "smearingWidth") {
        smearingWidth = parseDoubleWithUnits(val);
      }

      if (parameterName == "constantRelaxationTime") {
        constantRelaxationTime = parseDoubleWithUnits(val);
      }

      if (parameterName == "scatteringMatrixInMemory") {
        scatteringMatrixInMemory = parseBool(val);
      }

      if (parameterName == "useSymmetries") {
        useSymmetries = parseBool(val);
      }

      if (parameterName == "withIsotopeScattering") {
        withIsotopeScattering = parseBool(val);
      }

      if (parameterName == "massVariance") {
        std::vector<double> x = parseDoubleList(val);
        massVariance = Eigen::VectorXd::Zero(x.size());
        for (long unsigned i = 0; i < x.size(); i++) {
          massVariance(i) = x[i];
        }
      }

      if (parameterName == "boundaryLength") {
        boundaryLength = parseDoubleWithUnits(val);
      }

      // EPA
      if (parameterName == "epaFileName") {
        epaFileName = parseString(val);
      }

      if (parameterName == "minChemicalPotential") {
        minChemicalPotential = parseDoubleWithUnits(val);
      }

      if (parameterName == "maxChemicalPotential") {
        maxChemicalPotential = parseDoubleWithUnits(val);
      }

      if (parameterName == "deltaChemicalPotential") {
        deltaChemicalPotential = parseDoubleWithUnits(val);
      }

      if (parameterName == "minTemperature") {
        minTemperature = parseDouble(val) / temperatureAuToSi;
      }

      if (parameterName == "maxTemperature") {
        maxTemperature = parseDouble(val) / temperatureAuToSi;
      }

      if (parameterName == "deltaTemperature") {
        deltaTemperature = parseDouble(val) / temperatureAuToSi;
      }

      if (parameterName == "eFermiRange") {
        eFermiRange = parseDoubleWithUnits(val);
      }

      if (parameterName == "epaSmearingEnergy") {
        epaSmearingEnergy = parseDoubleWithUnits(val);
      }
      if (parameterName == "epaDeltaEnergy") {
        epaDeltaEnergy = parseDoubleWithUnits(val);
      }
      if (parameterName == "epaNumBins") {
        epaNumBins = parseLong(val);
      }
      if (parameterName == "epaMinEnergy") {
        epaMinEnergy = parseDoubleWithUnits(val);
      }
      if (parameterName == "epaMaxEnergy") {
        epaMaxEnergy = parseDoubleWithUnits(val);
      }
      if (parameterName == "epaEnergyRange") {
        epaEnergyRange = parseDoubleWithUnits(val);
      }
      if (parameterName == "epaEnergyStep") {
        epaEnergyStep = parseDoubleWithUnits(val);
      }

      // ELPH coupling plot App

      if (parameterName == "g2PlotStyle") {
        g2PlotStyle = parseString(val);
      }

      if (parameterName == "g2FixedPoint") {
        std::vector<double> x = parseDoubleList(val);
        for (auto i : {0, 1, 2}) {
          g2PlotFixedPoint(i) = x[i];
        }
      }

      if (parameterName == "g2PlotBandEl1") {
        std::vector<int> x = parseIntList(val);
        g2PlotEl1Bands.first = x[0];
        g2PlotEl1Bands.second = x[1];
      }

      if (parameterName == "g2PlotBandEl2") {
        std::vector<int> x = parseIntList(val);
        g2PlotEl2Bands.first = x[0];
        g2PlotEl2Bands.second = x[1];
      }

      if (parameterName == "g2PlotBandPh") {
        std::vector<int> x = parseIntList(val);
        g2PlotPhBands.first = x[0];
        g2PlotPhBands.second = x[1];
      }

      //////////////////////////////////////////

    } else { // it might be a block, or its content

      auto tup = parseBlockNameValue(lines, lineCounter);
      auto blockName = std::get<0>(tup);
      auto value = std::get<1>(tup);

      if (blockName == "crystal") {
        auto tup1 = parseCrystal(value);
        auto inputAtomicPositions_ = std::get<0>(tup1);
        auto inputAtomicSpecies_ = std::get<1>(tup1);
        auto inputSpeciesNames_ = std::get<2>(tup1);
        inputAtomicPositions = inputAtomicPositions_;
        inputAtomicSpecies = inputAtomicSpecies_;
        inputSpeciesNames = inputSpeciesNames_;
      }
      if (blockName == "point path") {
        auto tup = parsePathExtrema(value);
        pathLabels = std::get<0>(tup);
        pathExtrema = std::get<1>(tup);
      }
    }

    ///////////////////////////////////////////////

    lineCounter += 1;
  }
};

std::string Context::getPhD2FileName() { return phD2FileName; }
void Context::setPhD2FileName(const std::string x) { phD2FileName = x; }

std::string Context::getPhD3FileName() { return phD3FileName; }
void Context::setPhD3FileName(const std::string x) { phD3FileName = x; }

std::string Context::getSumRuleD2() { return sumRuleD2; }
void Context::setSumRuleD2(const std::string x) { sumRuleD2 = x; }

std::string Context::getEpwFileName() { return epwFileName; }
void Context::setEpwFileName(const std::string x) { epwFileName = x; }

std::string Context::getElectronH0Name() { return electronH0Name; }

void Context::setElectronH0Name(const std::string x) { electronH0Name = x; }

std::string Context::getWannier90Prefix() { return wannier90Prefix; }
void Context::setWannier90Prefix(const std::string x) { wannier90Prefix = x; }
std::string Context::getQuantumEspressoPrefix() {
  return quantumEspressoPrefix;
}
void Context::setQuantumEspressoPrefix(const std::string x) {
  quantumEspressoPrefix = x;
}

std::string Context::getElPhInterpolation() { return elPhInterpolation; }

double Context::getEpaSmearingEnergy() { return epaSmearingEnergy; }
double Context::getEpaDeltaEnergy() { return epaDeltaEnergy; }
double Context::getEpaMinEnergy() { return epaMinEnergy; }
double Context::getEpaMaxEnergy() { return epaMaxEnergy; }
int Context::getEpaNumBins() { return epaNumBins; }
double Context::getEpaEnergyRange() {return epaEnergyRange;}
double Context::getEpaEnergyStep() {return epaEnergyStep;}

double Context::getElectronFourierCutoff() { return electronFourierCutoff; }

std::string Context::getAppName() { return appName; }

Eigen::Vector3i Context::getQMesh() { return qMesh; }

Eigen::Vector3i Context::getKMesh() { return kMesh; }

std::string Context::getWindowType() { return windowType; }
void Context::setWindowType(const std::string x) { windowType = x; }

Eigen::Vector2d Context::getWindowEnergyLimit() { return windowEnergyLimit; }

void Context::setWindowEnergyLimit(const Eigen::Vector2d x) {
  windowEnergyLimit = x;
}

double Context::getWindowPopulationLimit() { return windowPopulationLimit; }
void Context::setWindowPopulationLimit(const double x) {
  windowPopulationLimit = x;
}

Eigen::VectorXd Context::getChemicalPotentials() { return chemicalPotentials; }

Eigen::VectorXd Context::getDopings() { return dopings; }

void Context::setDopings(const Eigen::VectorXd x) { dopings = x; }

Eigen::VectorXd Context::getTemperatures() { return temperatures; }

void Context::setTemperatures(const Eigen::VectorXd x) { temperatures = x; }

std::vector<std::string> Context::getSolverBTE() { return solverBTE; }

double Context::getConvergenceThresholdBTE() { return convergenceThresholdBTE; }

long Context::getMaxIterationsBTE() { return maxIterationsBTE; }

long Context::getDimensionality() { return dimensionality; }

double Context::getDosMinEnergy() { return dosMinEnergy; }

double Context::getDosMaxEnergy() { return dosMaxEnergy; }

double Context::getDosDeltaEnergy() { return dosDeltaEnergy; }

Eigen::MatrixXd Context::getInputAtomicPositions() {
  return inputAtomicPositions;
}

Eigen::VectorXi Context::getInputAtomicSpecies() { return inputAtomicSpecies; }

std::vector<std::string> Context::getInputSpeciesNames() {
  return inputSpeciesNames;
}

void Context::setInputAtomicPositions(const Eigen::MatrixXd x) {
  inputAtomicPositions = x;
}
void Context::setInputAtomicSpecies(const Eigen::VectorXi x) {
  inputAtomicSpecies = x;
}
void Context::setInputSpeciesNames(const std::vector<std::string> x) {
  inputSpeciesNames = x;
}

Eigen::Tensor<double, 3> Context::getPathExtrema() { return pathExtrema; }
std::vector<std::string> Context::getPathLabels() { return pathLabels; }

double Context::getDeltaPath() { return deltaPath; }

double Context::getFermiLevel() { return fermiLevel; }

void Context::setFermiLevel(const double &x) { fermiLevel = x; }

double Context::getNumOccupiedStates() { return numOccupiedStates; }

void Context::setNumOccupiedStates(const double &x) { numOccupiedStates = x; }

bool Context::getHasSpinOrbit() { return hasSpinOrbit; }

void Context::setHasSpinOrbit(const bool &x) { hasSpinOrbit = x; }

int Context::getSmearingMethod() { return smearingMethod; }

double Context::getSmearingWidth() { return smearingWidth; }
void Context::setSmearingWidth(const double x) { smearingWidth = x; }

double Context::getConstantRelaxationTime() { return constantRelaxationTime; }

bool Context::getScatteringMatrixInMemory() { return scatteringMatrixInMemory; }
void Context::setScatteringMatrixInMemory(const bool &x) {
  scatteringMatrixInMemory = x;
}

bool Context::getUseSymmetries() { return useSymmetries; }
void Context::setUseSymmetries(const bool &x) {
  useSymmetries = x;
}

Eigen::VectorXd Context::getMassVariance() { return massVariance; }

bool Context::getWithIsotopeScattering() { return withIsotopeScattering; }

double Context::getBoundaryLength() { return boundaryLength; }

std::string Context::getEpaFileName() {return epaFileName;}

double Context::getMinChemicalPotential() {return minChemicalPotential;}

double Context::getMaxChemicalPotential() {return maxChemicalPotential;}

double Context::getDeltaChemicalPotential() {return deltaChemicalPotential;}

double Context::getMinTemperature() {return minTemperature;}

double Context::getMaxTemperature() {return maxTemperature;}

double Context::getDeltaTemperature() {return deltaTemperature;}

double Context::getEFermiRange() {return eFermiRange;}

std::string Context::getG2PlotStyle() { return g2PlotStyle; }
void Context::setG2PlotStyle(const std::string x) { g2PlotStyle = x; }

Eigen::Vector3d Context::getG2PlotFixedPoint() { return g2PlotFixedPoint; }
void Context::setG2PlotFixedPoint(const Eigen::Vector3d x) {
  g2PlotFixedPoint = x;
}

std::pair<int, int> Context::getG2PlotEl1Bands() { return g2PlotEl1Bands; }
void Context::setG2PlotEl1Bands(const std::pair<int, int> x) {
  g2PlotEl1Bands = x;
}

std::pair<int, int> Context::getG2PlotEl2Bands() { return g2PlotEl2Bands; }
void Context::setG2PlotEl2Bands(const std::pair<int, int> x) {
  g2PlotEl2Bands = x;
}

std::pair<int, int> Context::getG2PlotPhBands() { return g2PlotPhBands; }
void Context::setG2PlotPhBands(const std::pair<int, int> x) {
  g2PlotPhBands = x;
}

