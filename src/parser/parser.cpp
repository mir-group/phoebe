#include "parser.h"

std::tuple<Crystal, PhononH0> Parser::parsePhHarmonic(Context &context) {

  std::string fileName = context.getPhonopyDispFileName();

  // check if this file is set -- if it is, we read from phonopy
  if (fileName.empty()) {
    if(context.getPhFC2FileName().find("hdf5") != std::string::npos) {
      // Note for the when someone forgets to set PhonopyDispFileName and gets an stoi error
      if(mpi->mpiHead()) {
        std::cout << "Warning: Reading in a QE .fc file with a name containing hdf5." << '\n' <<
        "If you are trying to use Phonopy force constants," <<
        " make sure to set the PhonopyDispFileName variable." << std::endl;
        }
    }
    return QEParser::parsePhHarmonic(context);
  }
  else {
    return PhonopyParser::parsePhHarmonic(context);
  }
}

std::tuple<Crystal, ElectronH0Fourier> Parser::parseElHarmonicFourier(
        Context &context) {
  return QEParser::parseElHarmonicFourier(context);
}

std::tuple<Crystal, ElectronH0Wannier> Parser::parseElHarmonicWannier(
            Context &context, Crystal *inCrystal) {
  return QEParser::parseElHarmonicWannier(context, inCrystal);
}
