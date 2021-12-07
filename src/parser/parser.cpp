#include "parser.h"

std::tuple<Crystal, PhononH0> Parser::parsePhHarmonic(Context &context) {

  std::string fileName = context.getPhonopyDispFileName();

  // check if this file is set -- if it is, we read from phonopy
  if (fileName.empty()) {
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
